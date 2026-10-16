# CCkit 功能说明：
"""Interactive Windows visualization shell for CCkit feature scripts.

The viewer uses PySide6 and PyVista/VTK. It reads a crystal structure and
``gbin.json``, samples a CNID grid or scans CSL candidates in worker threads,
and renders the selected frame with a hardware-accelerated VTK viewport.
Extended XYZ and ASE trajectories can also be opened and exported.
Feature names and their script ownership are supplied by ``feature_registry``;
OmniGB is the first connected card and future scripts can add cards without
replacing the application shell. The application is intentionally separate
from the CCkit command dispatcher: it is a desktop inspection tool rather
than a batch command.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import csv
import faulthandler
import importlib.util
import io
import json
import math
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator

from feature_registry import ViewerFeature, create_default_registry


if getattr(sys, "frozen", False):
    # PyInstaller stores bundled data below _MEIPASS.  Keeping this lookup in
    # one place lets the same source run both from the repository and as an
    # onedir executable.
    REPOSITORY_ROOT = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
else:
    REPOSITORY_ROOT = Path(__file__).resolve().parent
CORE_DIRECTORY = REPOSITORY_ROOT / "core"
FEATURE_REGISTRY = create_default_registry(REPOSITORY_ROOT, core_directory=CORE_DIRECTORY)
MODEL_SCRIPT = FEATURE_REGISTRY.get("omnigb").script_path
STRUCTURE_FILE_FILTER = (
    "Structure files ("
    "*.cif *.vasp POSCAR* CONTCAR* *.xyz *.pdb *.res *.cssr *.xsf "
    "*.mson *.json *.yaml *.yml);;"
    "All files (*)"
)
STARTUP_LOG = Path(tempfile.gettempdir()) / "CCkitGBViewer_startup.log"


@dataclass(frozen=True)
class StructureExportRequest:
    """Keep the chosen file and model settings together during an export."""

    target: Path
    fmt: str
    input_path: Path
    config: dict[str, Any]
    feature: ViewerFeature
    scope: str = "current"
    frame_index: int = 0


@dataclass(frozen=True)
class StructureTrajectory:
    """Keep independently built frames and their RBT coordinates together."""

    frames: tuple[Any, ...]
    metadata: tuple[dict[str, Any], ...]
    details: dict[str, Any] | None = None


def _viewer_config(config: dict[str, Any]) -> dict[str, Any]:
    """Use CNID sampling in the GUI; the CLI retains legacy translations."""

    result = copy.deepcopy(config)
    result["cnid_grid"] = list(_load_model_module().cnid_grid(result))
    if any(value > 1000 for value in result["cnid_grid"]):
        raise ValueError("the viewer supports 1–1000 subdivisions per CNID direction")
    result.pop("translation", None)
    result.pop("dydz", None)
    return result


def _read_trajectory(path: Path) -> StructureTrajectory:
    """Load all lattice frames from extended XYZ or ASE trajectory files."""

    import numpy as np
    try:
        from ase.io import iread
    except ImportError as exc:
        raise RuntimeError(
            "Trajectory I/O requires ASE. Install requirements-gb-viewer.txt "
            "in the Python environment running this viewer."
        ) from exc
    from pymatgen.core import Lattice, Structure

    frames = []
    metadata = []
    for index, atoms in enumerate(iread(str(path), index=":")):
        if len(atoms) == 0:
            raise ValueError(f"frame {index} contains no atoms")
        if not np.isfinite(atoms.cell.array).all() or abs(atoms.cell.volume) < 1e-12:
            raise ValueError(
                f"frame {index} has no valid lattice. Use extended XYZ with "
                "a Lattice field, or an ASE .traj file."
            )
        lattice = Lattice(atoms.cell.array, pbc=tuple(bool(x) for x in atoms.pbc))
        frames.append(Structure(
            lattice, atoms.get_chemical_symbols(), atoms.positions,
            coords_are_cartesian=True,
        ))
        metadata.append({**atoms.info, "frame": index})
    if not frames:
        raise ValueError("the trajectory contains no frames")
    return StructureTrajectory(tuple(frames), tuple(metadata))


def _write_xyz_frames(
    trajectory: StructureTrajectory, indices: range | list[int], target: Path
) -> None:
    """Write extended XYZ with the cell, PBC and RBT sampling coordinates."""

    try:
        from ase import Atoms
        from ase.io import write
    except ImportError as exc:
        raise RuntimeError(
            "XYZ export requires ASE. Install requirements-gb-viewer.txt "
            "in the Python environment running this viewer."
        ) from exc

    def images() -> Iterator[Any]:
        for index in indices:
            structure = trajectory.frames[index]
            atoms = Atoms(
                symbols=[str(site.specie) for site in structure],
                positions=structure.cart_coords,
                cell=structure.lattice.matrix,
                pbc=structure.lattice.pbc,
            )
            atoms.info.update(trajectory.metadata[index])
            atoms.info["frame"] = index
            yield atoms

    write(str(target), images(), format="extxyz", write_results=False)


def _export_trajectory(
    trajectory: StructureTrajectory, request: StructureExportRequest
) -> None:
    """Write a saved export snapshot and replace the destination on success."""

    indices: range | list[int]
    if request.scope == "all":
        if request.fmt != "extxyz":
            raise ValueError("all frames must be exported as extended XYZ")
        indices = range(len(trajectory.frames))
    else:
        if not 0 <= request.frame_index < len(trajectory.frames):
            raise ValueError("the selected frame is outside the current CNID grid")
        indices = [request.frame_index]
    request.target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".cckit_export_", dir=request.target.parent
    ) as directory:
        temporary = Path(directory) / request.target.name
        if request.fmt == "extxyz":
            _write_xyz_frames(trajectory, indices, temporary)
        else:
            trajectory.frames[request.frame_index].to(
                fmt=request.fmt, filename=str(temporary)
            )
        temporary.replace(request.target)


def _configure_diagnostics() -> Any:
    """Keep Python and native crash details for double-click launches."""

    stream = STARTUP_LOG.open("w", encoding="utf-8", buffering=1)
    stream.write(f"Starting CCkit GB Viewer: {sys.executable}\n")
    if sys.stdout is None:
        sys.stdout = stream
    if sys.stderr is None:
        sys.stderr = stream
    faulthandler.enable(file=stream)
    return stream


def _load_model_module(script_path: Path = MODEL_SCRIPT) -> Any:
    """Load the existing model-gb implementation without copying its logic."""

    if not script_path.is_file():
        raise FileNotFoundError(f"feature script was not found: {script_path}")
    spec = importlib.util.spec_from_file_location("cckit_feature_script", script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"could not load feature script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_config(path: Path) -> dict[str, Any]:
    """Read a JSON object used by model-gb."""

    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("gbin.json must contain a JSON object")
    if data.get("mode", "build") not in {"build", "scan"}:
        raise ValueError("mode must be either 'build' or 'scan'")
    return _viewer_config(data)


def _write_config(path: Path, config: dict[str, Any]) -> None:
    """Write a formatted gbin.json while preserving Unicode text."""

    path.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _vector(config: dict[str, Any], key: str, default: list[float]) -> list[float]:
    """Return a numeric three-component vector for a viewer control."""

    value = config.get(key, default)
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{key} must be a list")
    if len(value) != 3:
        raise ValueError(f"{key} must contain three numbers")
    try:
        result = [float(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must contain only numbers") from exc
    if not all(math.isfinite(item) for item in result):
        raise ValueError(f"{key} must contain only finite numbers")
    return result


def _termination_values(config: dict[str, Any]) -> tuple[float, float]:
    """Read the same nested or flat dp settings accepted by model-gb."""

    shifts = config.get("termination_shift", {})
    if shifts is None:
        shifts = {}
    if not isinstance(shifts, dict):
        raise ValueError("termination_shift must be a JSON object")
    return (
        float(config.get("dp1", shifts.get("dp1", 0.0))),
        float(config.get("dp2", shifts.get("dp2", 0.0))),
    )


def _build_structure(
    input_path: Path,
    config: dict[str, Any],
    feature: ViewerFeature | None = None,
    *,
    progress: Callable[[int, int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> StructureTrajectory:
    """Build preview frames using the registered feature's own backend."""

    feature = feature or FEATURE_REGISTRY.get("omnigb")
    if not feature.available or feature.script_path is None:
        raise ValueError(f"feature is not implemented yet: {feature.name}")
    model_gb = _load_model_module(feature.script_path)
    with contextlib.redirect_stdout(io.StringIO()):
        if feature.key == "wyckoff":
            builder = model_gb.build_wyckoff_frames
        elif feature.key == "interface":
            builder = model_gb.build_interface_frames
        else:
            builder = model_gb.build_boundary_frames
        result = builder(input_path, config, progress=progress, cancelled=cancelled)
    return StructureTrajectory(
        tuple(result["frames"]), tuple(result["metadata"]), result.get("report")
    )


def _scan_structure(
    input_path: Path,
    config: dict[str, Any],
    feature: ViewerFeature,
) -> tuple[list[dict[str, str]], str]:
    """Run the existing CSL scan and keep its successful rows and CSV text."""

    if not feature.available or feature.script_path is None:
        raise ValueError(f"feature is not implemented yet: {feature.name}")
    model_gb = _load_model_module(feature.script_path)
    with tempfile.TemporaryDirectory(prefix="cckit_gb_scan_") as work_dir:
        output_path = Path(work_dir) / "GB_scan_results.csv"
        with contextlib.redirect_stdout(io.StringIO()):
            model_gb.scan_candidates(input_path, config, output_path)
        csv_text = output_path.read_text(encoding="utf-8")
        rows = list(csv.DictReader(io.StringIO(csv_text)))
    return rows, csv_text


SCAN_COLUMNS = (
    "angle_deg",
    "sigma_1",
    "sigma_2",
    "deformation_max_abs",
    "csl_mode",
    "status",
)


def _format_deformation(value: str) -> str:
    """Format deformation without hiding values that would round to zero."""

    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number != 0.0 and abs(number) < 0.01:
        return f"{number:.2e}"
    return f"{number:.2f}"


def _format_scan_table(rows: list[dict[str, str]]) -> str:
    """Return a fixed-width view of successful scan rows for the GUI."""

    formatted_rows = []
    for row in rows:
        formatted_rows.append(
            [
                _format_deformation(row.get(column, ""))
                if column == "deformation_max_abs"
                else str(row.get(column, ""))
                for column in SCAN_COLUMNS
            ]
        )
    header = list(SCAN_COLUMNS)
    widths = [len(column) for column in SCAN_COLUMNS]
    for row in formatted_rows:
        widths = [max(width, len(value)) for width, value in zip(widths, row)]
    lines = [
        "  ".join(column.ljust(width) for column, width in zip(header, widths)),
        "  ".join("-" * width for width in widths),
    ]
    lines.extend(
        "  ".join(value.ljust(width) for value, width in zip(row, widths))
        for row in formatted_rows
    )
    return "\n".join(lines)


def _omnigb_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Read UI values with the same aliases and defaults as model-gb."""

    model_gb = _load_model_module()
    get_value = model_gb._config_value
    mode = config.get("mode", "build")
    if mode not in {"build", "scan"}:
        raise ValueError("mode must be either 'build' or 'scan'")
    axis_input = get_value(config, ("rotation_axis", "axis"), [0, 0, 1])
    axis_notation = model_gb.index_notation(axis_input, config.get("rotation_axis_notation"))
    axis = model_gb.index_values(axis_input, "direction", axis_notation)
    plane_input = get_value(config, ("gb_plane", "grain_boundary_plane", "plane", "hkl"), [1, 0, 0])
    plane_notation = model_gb.index_notation(plane_input, config.get("gb_plane_notation"))
    plane = model_gb.index_values(plane_input, "plane", plane_notation)
    if mode == "build" and not any(plane):
        raise ValueError("gb_plane must not be the zero vector")
    sigma = model_gb._positive_int(
        get_value(config, ("sigma", "target_sigma"), 5), "sigma"
    )
    angle_range = get_value(
        config, ("angle_range_deg", "angle_range", "theta_range"), [0.0, 180.0]
    )
    if not isinstance(angle_range, (list, tuple)) or len(angle_range) != 2:
        raise ValueError("angle_range_deg must contain [start_deg, stop_deg]")
    values = {
        "mode": mode,
        "rotation_axis": axis,
        "rotation_axis_notation": axis_notation,
        "gb_plane": plane,
        "gb_plane_notation": plane_notation,
        "sigma": sigma,
        "max_sigma": model_gb._positive_int(
            get_value(config, ("max_sigma", "maxsigma"), sigma), "max_sigma"
        ),
        "rotation_angle_deg": model_gb._finite_float(
            get_value(
                config,
                ("rotation_angle_deg", "angle_deg", "rotation_angle", "theta"),
                0.0,
            ),
            "rotation_angle_deg",
        ),
        "angle_start_deg": model_gb._finite_float(
            config.get("angle_start_deg", angle_range[0]), "angle_start_deg"
        ),
        "angle_stop_deg": model_gb._finite_float(
            config.get("angle_stop_deg", angle_range[1]), "angle_stop_deg"
        ),
        "angle_step_deg": model_gb._finite_float(
            get_value(config, ("angle_step_deg", "angle_step"), 1.0),
            "angle_step_deg",
        ),
        "csl_mode": model_gb._csl_mode(config)[0],
        "tolerances": model_gb._tolerances(config),
    }
    if values["angle_step_deg"] <= 0:
        raise ValueError("angle_step_deg must be greater than zero")
    if values["angle_stop_deg"] < values["angle_start_deg"]:
        raise ValueError(
            "angle_stop_deg must be greater than or equal to angle_start_deg"
        )
    return values


def _qt_imports() -> tuple[Any, ...]:
    """Import optional GUI dependencies and provide one actionable error."""

    try:
        import numpy as np
        import pyvista as pv
        from pyvistaqt import QtInteractor
        from PySide6.QtCore import QSignalBlocker, QThread, QTimer, Qt, Signal
        from PySide6.QtGui import QAction
        from PySide6.QtWidgets import (
            QApplication,
            QComboBox,
            QFileDialog,
            QFormLayout,
            QDoubleSpinBox,
            QGridLayout,
            QGroupBox,
            QHBoxLayout,
            QLabel,
            QMainWindow,
            QPlainTextEdit,
            QPushButton,
            QScrollArea,
            QSizePolicy,
            QSplitter,
            QStackedWidget,
            QVBoxLayout,
            QWidget,
        )
    except ImportError as exc:
        raise RuntimeError(
            "GB viewer requires PySide6, pyvista, pyvistaqt, and vtk. "
            "Install them with: python -m pip install -r "
            "requirements-gb-viewer.txt "
            f"(import error: {exc})"
        ) from exc

    return (
        np,
        pv,
        QtInteractor,
        QSignalBlocker,
        QThread,
        QTimer,
        Qt,
        Signal,
        QAction,
        QApplication,
        QComboBox,
        QFileDialog,
        QFormLayout,
        QDoubleSpinBox,
        QGridLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QPlainTextEdit,
        QPushButton,
        QScrollArea,
        QSizePolicy,
        QSplitter,
        QStackedWidget,
        QVBoxLayout,
        QWidget,
    )


def _report_startup_error(error: BaseException) -> None:
    """Write a startup traceback and show it when the GUI has no console."""

    message = "CCkit GB Viewer could not start.\n\n" + "".join(
        traceback.format_exception(type(error), error, error.__traceback__)
    )
    try:
        with STARTUP_LOG.open("a", encoding="utf-8") as stream:
            stream.write(message)
    except OSError:
        pass
    try:
        import ctypes

        ctypes.windll.user32.MessageBoxW(
            0,
            f"{message}\n\nLog: {STARTUP_LOG}",
            "CCkit Grain Boundary Viewer",
            0x10,
        )
    except (AttributeError, OSError):
        print(message, file=sys.stderr)


def _lattice_lines(np: Any, pv: Any, lattice: Any) -> Any:
    """Return a PyVista line mesh for the eight unit-cell corners."""

    a, b, c = np.asarray(lattice.matrix, dtype=float)
    corners = np.array(
        [
            [0, 0, 0],
            a,
            b,
            c,
            a + b,
            a + c,
            b + c,
            a + b + c,
        ],
        dtype=float,
    )
    pairs = (
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (1, 5),
        (2, 4),
        (2, 6),
        (3, 5),
        (3, 6),
        (4, 7),
        (5, 7),
        (6, 7),
    )
    lines = np.asarray([[2, start, end] for start, end in pairs], dtype=int)
    return pv.PolyData(corners, lines=lines.ravel())


def _bond_lines(
    np: Any,
    pv: Any,
    structure: Any,
    minimum: float,
    maximum: float,
    pair_rules: dict[str, dict[str, Any]] | None = None,
) -> Any:
    """Create periodic bond segments, optionally using element-pair rules."""

    search_maximum = maximum
    for rule in (pair_rules or {}).values():
        try:
            search_maximum = max(search_maximum, float(rule.get("maximum", maximum)))
        except (AttributeError, TypeError, ValueError):
            continue
    center, neighbour, images, distances = structure.get_neighbor_list(search_maximum)
    points: list[Any] = []
    line_cells: list[int] = []
    seen: set[tuple[int, int, tuple[int, int, int]]] = set()
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    coordinates = np.asarray(structure.cart_coords, dtype=float)
    species = [str(site.specie) for site in structure]
    for first, second, image, distance in zip(center, neighbour, images, distances):
        pair_name = "-".join(sorted((species[int(first)], species[int(second)])))
        rule = (pair_rules or {}).get(pair_name)
        if rule is not None:
            if not bool(rule.get("enabled", True)):
                continue
            pair_minimum = float(rule.get("minimum", minimum))
            pair_maximum = float(rule.get("maximum", maximum))
        else:
            pair_minimum, pair_maximum = minimum, maximum
        if float(distance) < pair_minimum or float(distance) > pair_maximum:
            continue
        first = int(first)
        second = int(second)
        image_tuple = tuple(int(value) for value in image)
        if first == second:
            if image_tuple == (0, 0, 0):
                continue
            canonical_image = min(image_tuple, tuple(-value for value in image_tuple))
            key = (first, second, canonical_image)
        elif first < second:
            key = (first, second, image_tuple)
        else:
            key = (second, first, tuple(-value for value in image_tuple))
        if key in seen:
            continue
        seen.add(key)
        endpoint_a = coordinates[first]
        endpoint_b = coordinates[second] + np.dot(image_tuple, lattice)
        points.extend((endpoint_a, endpoint_b))
        start = len(points) - 2
        line_cells.extend((2, start, start + 1))
    if not points:
        return pv.PolyData()
    return pv.PolyData(np.asarray(points), lines=np.asarray(line_cells, dtype=int))


def _grain_boundary_planes(
    np: Any, pv: Any, structure: Any, fractions: Any = (0.0, 0.5)
) -> tuple[Any, ...]:
    """Create the two periodic interface planes in the final GB cell."""

    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    meshes = []
    for normal_fraction in fractions:
        fractional = np.asarray(
            [
                [normal_fraction, 0.0, 0.0],
                [normal_fraction, 1.0, 0.0],
                [normal_fraction, 1.0, 1.0],
                [normal_fraction, 0.0, 1.0],
            ],
            dtype=float,
        )
        vertices = np.dot(fractional, lattice)
        meshes.append(pv.PolyData(vertices, faces=np.asarray([4, 0, 1, 2, 3])))
    return tuple(meshes)


def create_application_classes() -> tuple[Any, Any]:
    """Create Qt classes after optional dependencies have been imported."""

    (
        np,
        pv,
        QtInteractor,
        QSignalBlocker,
        QThread,
        QTimer,
        Qt,
        Signal,
        QAction,
        QApplication,
        QComboBox,
        QFileDialog,
        QFormLayout,
        QDoubleSpinBox,
        QGridLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QPlainTextEdit,
        QPushButton,
        QScrollArea,
        QSizePolicy,
        QSplitter,
        QStackedWidget,
        QVBoxLayout,
        QWidget,
    ) = _qt_imports()

    from PySide6.QtCore import QSettings
    from PySide6.QtGui import QColor
    from PySide6.QtWidgets import (
        QCheckBox,
        QColorDialog,
        QDialog,
        QDialogButtonBox,
        QToolBar,
        QSlider,
        QSpinBox,
    )

    class NoWheelDoubleSpinBox(QDoubleSpinBox):
        """Consume wheel events so numeric controls cannot react to them."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)

        def wheelEvent(self, event: Any) -> None:
            event.accept()

    class NoWheelSpinBox(QSpinBox):
        """Keep integer editors from changing values during panel scrolling."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)

        def wheelEvent(self, event: Any) -> None:
            event.accept()

    class NoWheelSlider(QSlider):
        """Change frames by dragging, clicking or keys, without wheel events."""

        def wheelEvent(self, event: Any) -> None:
            event.accept()

    class NoWheelComboBox(QComboBox):
        """Consume wheel events so option controls cannot change selection."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)

        def wheelEvent(self, event: Any) -> None:
            event.accept()

    # All numeric and option editors in the viewer use the same non-wheel-
    # interactive behavior, including personal preferences.
    QDoubleSpinBox = NoWheelDoubleSpinBox
    QComboBox = NoWheelComboBox
    QSpinBox = NoWheelSpinBox
    QSlider = NoWheelSlider

    class FeatureCardStack(QWidget):
        """Keep independent feature pages in a top-aligned, visible-only layout."""

        def __init__(self) -> None:
            super().__init__()
            self._pages: list[Any] = []
            self._current_index = -1
            self._page_layout = QVBoxLayout(self)
            self._page_layout.setContentsMargins(0, 0, 0, 0)
            self._page_layout.setSpacing(0)
            self._page_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
            self.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)

        def addWidget(self, page: Any) -> int:
            """Retain a page's state while excluding it from layout when hidden."""

            index = len(self._pages)
            self._pages.append(page)
            self._page_layout.addWidget(page)
            page.hide()
            if self._current_index < 0:
                self.setCurrentIndex(index)
            return index

        def setCurrentIndex(self, index: int) -> None:
            """Show only the selected page directly below the feature selector."""

            if not 0 <= index < len(self._pages):
                raise IndexError(f"feature page index out of range: {index}")
            for page_index, page in enumerate(self._pages):
                page.setVisible(page_index == index)
            self._current_index = index
            # A normal box layout excludes explicitly hidden pages, including
            # their minimum height and height-for-width calculations.
            self._page_layout.invalidate()
            self.updateGeometry()

    class CrystalIndexEditor(QWidget):
        """Switch three/four indices while preserving a physical axis or plane."""

        changed = Signal()
        invalid = Signal(str)

        def __init__(self, kind: str) -> None:
            super().__init__()
            self.kind = kind
            self._mode = "three"
            self._updating = True
            self._backend = _load_model_module()
            row = QHBoxLayout(self)
            row.setContentsMargins(0, 0, 0, 0)
            self.labels = []
            self.spins = []
            for _ in range(4):
                label = QLabel()
                spin = QSpinBox()
                spin.setRange(-2000000, 2000000)
                spin.setMaximumWidth(90)
                spin.setKeyboardTracking(False)
                spin.valueChanged.connect(self._values_changed)
                row.addWidget(label)
                row.addWidget(spin)
                self.labels.append(label)
                self.spins.append(spin)
            self.selector = QComboBox()
            self.selector.addItem("三指数", "three")
            self.selector.addItem("六方四指数", "four")
            self.selector.currentIndexChanged.connect(self._notation_changed)
            row.addWidget(self.selector)
            self.set_indices([0, 0, 1] if kind == "direction" else [1, 0, 0], "three")
            self._updating = False

        def notation(self) -> str:
            return self._mode

        def indices(self) -> list[int]:
            order = (0, 1, 2, 3) if self._mode == "four" else (0, 1, 3)
            return [self.spins[i].value() for i in order]

        def commit(self) -> None:
            for spin in self.spins:
                spin.interpretText()

        def set_indices(self, value: Any, notation: str | None = None) -> None:
            notation = self._backend.index_notation(value, notation)
            values = self._backend.index_values(value, self.kind, notation)
            if any(item != round(item) for item in values):
                raise ValueError("GUI index controls require integers; scale fractional indices first")
            supplied = values if notation == "four" else [values[0], values[1], 0, values[2]]
            if any(not spin.minimum() <= item <= spin.maximum()
                   for spin, item in zip(self.spins, supplied)):
                raise ValueError("index exceeds the GUI input range")
            self._updating = True
            try:
                self._mode = notation
                names = ("u", "v", "t", "w") if self.kind == "direction" else ("h", "k", "i", "l")
                if notation == "three":
                    names = ("U", "V", "", "W") if self.kind == "direction" else ("h", "k", "", "l")
                for i, (spin, item) in enumerate(zip(self.spins, supplied)):
                    if not spin.minimum() <= item <= spin.maximum():
                        raise ValueError("index exceeds the GUI input range")
                    spin.setValue(int(item))
                    self.labels[i].setText(names[i])
                self.labels[2].setVisible(notation == "four")
                self.spins[2].setVisible(notation == "four")
                self.spins[2].setReadOnly(True)
                self.spins[2].setEnabled(False)
                blocker = QSignalBlocker(self.selector)
                self.selector.setCurrentIndex(self.selector.findData(notation))
                del blocker
            finally:
                self._updating = False

        def _values_changed(self, *_args: Any) -> None:
            if self._updating:
                return
            if self._mode == "four":
                blocker = QSignalBlocker(self.spins[2])
                self.spins[2].setValue(-self.spins[0].value() - self.spins[1].value())
                del blocker
            self.changed.emit()

        def _notation_changed(self, *_args: Any) -> None:
            if self._updating:
                return
            requested = str(self.selector.currentData())
            try:
                values = self.indices()
                if requested == "four":
                    converted = self._backend.four_indices(values, self.kind)
                else:
                    converted = self._backend.three_indices(values, self.kind, "four")
                self.set_indices(converted, requested)
            except ValueError as exc:
                blocker = QSignalBlocker(self.selector)
                self.selector.setCurrentIndex(self.selector.findData(self._mode))
                del blocker
                self.invalid.emit(str(exc))
                return
            self.changed.emit()


    class PreferenceDialog(QDialog):
        """Edit persistent rendering preferences for the atom viewport."""

        def __init__(
            self,
            preferences: dict[str, Any],
            elements: list[str] | tuple[str, ...] = (),
            parent: Any = None,
        ) -> None:
            super().__init__(parent)
            self.setWindowTitle("Personal preferences")
            self.setModal(True)
            self.plane_color = QColor(preferences["plane_color"])
            layout = QVBoxLayout(self)
            form = QFormLayout()

            self.radius_scale = QDoubleSpinBox()
            self.radius_scale.setRange(0.1, 5.0)
            self.radius_scale.setDecimals(2)
            self.radius_scale.setSingleStep(0.1)
            self.radius_scale.setValue(preferences["atom_radius_scale"])
            form.addRow("Atom size scale", self.radius_scale)

            self.show_bonds = QCheckBox("Show bonds")
            self.show_bonds.setChecked(preferences["show_bonds"])
            form.addRow(self.show_bonds)

            self.bond_minimum = QDoubleSpinBox()
            self.bond_minimum.setRange(0.0, 20.0)
            self.bond_minimum.setDecimals(3)
            self.bond_minimum.setSingleStep(0.1)
            self.bond_minimum.setValue(preferences["bond_minimum"])
            form.addRow("Bond minimum (Å)", self.bond_minimum)

            self.bond_maximum = QDoubleSpinBox()
            self.bond_maximum.setRange(0.01, 20.0)
            self.bond_maximum.setDecimals(3)
            self.bond_maximum.setSingleStep(0.1)
            self.bond_maximum.setValue(preferences["bond_maximum"])
            form.addRow("Bond maximum (Å)", self.bond_maximum)

            self.pair_rows: dict[str, tuple[Any, Any, Any]] = {}
            pair_group = QGroupBox("Element-pair bond rules")
            pair_layout = QGridLayout(pair_group)
            pair_layout.addWidget(QLabel("Pair"), 0, 0)
            pair_layout.addWidget(QLabel("Enabled"), 0, 1)
            pair_layout.addWidget(QLabel("Minimum (Å)"), 0, 2)
            pair_layout.addWidget(QLabel("Maximum (Å)"), 0, 3)
            stored_rules = preferences.get("bond_pairs", {})
            if not isinstance(stored_rules, dict):
                stored_rules = {}
            unique_elements = sorted(set(elements))
            row = 1
            for first_index, first in enumerate(unique_elements):
                for second in unique_elements[first_index:]:
                    pair_name = "-".join(sorted((first, second)))
                    rule = stored_rules.get(pair_name, {})
                    if not isinstance(rule, dict):
                        rule = {}
                    pair_layout.addWidget(QLabel(pair_name), row, 0)
                    enabled = QCheckBox()
                    enabled.setChecked(bool(rule.get("enabled", True)))
                    pair_layout.addWidget(enabled, row, 1)
                    pair_minimum = QDoubleSpinBox()
                    pair_minimum.setRange(0.0, 20.0)
                    pair_minimum.setDecimals(3)
                    pair_minimum.setSingleStep(0.1)
                    pair_minimum.setValue(float(rule.get("minimum", preferences["bond_minimum"])))
                    pair_layout.addWidget(pair_minimum, row, 2)
                    pair_maximum = QDoubleSpinBox()
                    pair_maximum.setRange(0.01, 20.0)
                    pair_maximum.setDecimals(3)
                    pair_maximum.setSingleStep(0.1)
                    pair_maximum.setValue(float(rule.get("maximum", preferences["bond_maximum"])))
                    pair_layout.addWidget(pair_maximum, row, 3)
                    self.pair_rows[pair_name] = (enabled, pair_minimum, pair_maximum)
                    row += 1
            pair_hint = QLabel(
                "Each listed pair has its own range. Element pairs not listed "
                "use the global cutoff."
            )
            pair_hint.setWordWrap(True)
            pair_hint.setStyleSheet("color: #475569;")
            pair_layout.addWidget(pair_hint, row, 0, 1, 4)
            layout.addWidget(pair_group)

            self.show_plane = QCheckBox("Show symmetric GB planes")
            self.show_plane.setChecked(preferences["show_plane"])
            form.addRow(self.show_plane)

            self.plane_opacity = QDoubleSpinBox()
            self.plane_opacity.setRange(0.05, 1.0)
            self.plane_opacity.setDecimals(2)
            self.plane_opacity.setSingleStep(0.05)
            self.plane_opacity.setValue(preferences["plane_opacity"])
            form.addRow("Plane opacity", self.plane_opacity)

            self.plane_color_button = QPushButton()
            self.plane_color_button.clicked.connect(self.choose_plane_color)
            self._update_color_button()
            form.addRow("Plane color", self.plane_color_button)
            layout.addLayout(form)

            self.validation_label = QLabel()
            self.validation_label.setStyleSheet("color: #b42318;")
            layout.addWidget(self.validation_label)
            buttons = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok
                | QDialogButtonBox.StandardButton.Cancel
            )
            buttons.accepted.connect(self.accept_preferences)
            buttons.rejected.connect(self.reject)
            layout.addWidget(buttons)

        def _update_color_button(self) -> None:
            self.plane_color_button.setText(self.plane_color.name())
            self.plane_color_button.setStyleSheet(
                f"background-color: {self.plane_color.name()};"
            )

        def choose_plane_color(self) -> None:
            selected = QColorDialog.getColor(self.plane_color, self, "GB plane color")
            if selected.isValid():
                self.plane_color = selected
                self._update_color_button()

        def accept_preferences(self) -> None:
            if self.bond_minimum.value() >= self.bond_maximum.value():
                self.validation_label.setText(
                    "Bond maximum must be greater than bond minimum."
                )
                return
            for pair_name, (_, minimum, maximum) in self.pair_rows.items():
                if minimum.value() >= maximum.value():
                    self.validation_label.setText(
                        f"Bond maximum must be greater than minimum for {pair_name}."
                    )
                    return
            self.accept()

        def values(self) -> dict[str, Any]:
            """Return validated display preferences."""

            return {
                "atom_radius_scale": self.radius_scale.value(),
                "show_bonds": self.show_bonds.isChecked(),
                "bond_minimum": self.bond_minimum.value(),
                "bond_maximum": self.bond_maximum.value(),
                "bond_pairs": {
                    pair_name: {
                        "enabled": enabled.isChecked(),
                        "minimum": minimum.value(),
                        "maximum": maximum.value(),
                    }
                    for pair_name, (enabled, minimum, maximum) in self.pair_rows.items()
                },
                "show_plane": self.show_plane.isChecked(),
                "plane_opacity": self.plane_opacity.value(),
                "plane_color": self.plane_color.name(),
            }

    class StructureExportDialog(QDialog):
        """Choose one displayed frame or the complete sequence before saving."""

        def __init__(self, parent: Any, frame_index: int, frame_count: int) -> None:
            super().__init__(parent)
            self.setWindowTitle("Export structure / trajectory")
            layout = QFormLayout(self)
            self.scope_combo = QComboBox()
            self.scope_combo.addItem(f"Current frame ({frame_index})", "current")
            self.scope_combo.addItem(f"All frames ({frame_count})", "all")
            self.format_combo = QComboBox()
            self.scope_combo.currentIndexChanged.connect(self._update_formats)
            layout.addRow("Export", self.scope_combo)
            layout.addRow("File format", self.format_combo)
            self._update_formats()
            buttons = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            )
            buttons.accepted.connect(self.accept)
            buttons.rejected.connect(self.reject)
            layout.addRow(buttons)

        def _update_formats(self, *_args: Any) -> None:
            selected = self.format_combo.currentData()
            self.format_combo.clear()
            if self.scope_combo.currentData() == "current":
                self.format_combo.addItem("POSCAR", "poscar")
                self.format_combo.addItem("CIF", "cif")
            self.format_combo.addItem("Extended XYZ (*.xyz)", "extxyz")
            index = self.format_combo.findData(selected)
            if index >= 0:
                self.format_combo.setCurrentIndex(index)

    class BuildWorker(QThread):
        """Build a preview without blocking the Qt event loop."""

        completed = Signal(object)
        progress = Signal(int, int)
        failed = Signal(str)

        def __init__(
            self,
            input_path: Path,
            config: dict[str, Any],
            feature: ViewerFeature,
        ) -> None:
            super().__init__()
            self.input_path = input_path
            self.config = copy.deepcopy(config)
            self.feature = feature

        def run(self) -> None:
            try:
                structure = _build_structure(
                    self.input_path, self.config, self.feature,
                    progress=self.progress.emit, cancelled=self.isInterruptionRequested,
                )
            except Exception as exc:  # worker boundary must report all failures
                self.failed.emit(str(exc))
            else:
                self.completed.emit(structure)

    class StructureExportWorker(QThread):
        """Write a captured sequence without blocking frame navigation."""

        completed = Signal()
        failed = Signal(str)

        def __init__(
            self, trajectory: StructureTrajectory, request: StructureExportRequest,
            parent: Any,
        ) -> None:
            super().__init__(parent)
            self.trajectory = trajectory
            self.request = request

        def run(self) -> None:
            try:
                _export_trajectory(self.trajectory, self.request)
            except Exception as exc:
                self.failed.emit(str(exc))
            else:
                self.completed.emit()

    class ScanWorker(QThread):
        """Run a CSL scan without blocking the Qt event loop."""

        completed = Signal(object, str)
        failed = Signal(str)

        def __init__(
            self,
            input_path: Path,
            config: dict[str, Any],
            feature: ViewerFeature,
        ) -> None:
            super().__init__()
            self.input_path = input_path
            self.config = copy.deepcopy(config)
            self.feature = feature

        def run(self) -> None:
            try:
                rows, csv_text = _scan_structure(
                    self.input_path, self.config, self.feature
                )
            except Exception as exc:  # worker boundary must report all failures
                self.failed.emit(str(exc))
            else:
                self.completed.emit(rows, csv_text)

    class GrainBoundaryViewer(QMainWindow):
        """Desktop editor and renderer for model-gb build configurations."""

        def __init__(self, input_path: Path | None, config_path: Path | None) -> None:
            super().__init__()
            self.setWindowTitle("CCkit · Grain Boundary Viewer")
            self.resize(1440, 900)
            self.input_path = input_path
            self.config_path = config_path
            self.active_feature = FEATURE_REGISTRY.get("omnigb")
            self.worker: BuildWorker | None = None
            self.scan_worker: ScanWorker | None = None
            self.export_worker: StructureExportWorker | None = None
            self.last_structure: Any | None = None
            self.trajectory: StructureTrajectory | None = None
            self.current_frame_index = 0
            self.trajectory_source: Path | None = None
            self.rendered_config: dict[str, Any] | None = None
            self.rendered_input_path: Path | None = None
            self.rendered_feature_key: str | None = None
            self._export_request: StructureExportRequest | None = None
            self._export_build_active = False
            self._export_write_active = False
            self._closing = False
            self.scan_rows: list[dict[str, str]] = []
            self.scan_csv_text = ""
            self._camera_initialized = False
            self.preferences = self._load_preferences()
            self.build_serial = 0
            self.updating_controls = False
            self.auto_refresh = True
            self._build_timer = QTimer(self)
            self._build_timer.setSingleShot(True)
            self._build_timer.setInterval(350)
            self._build_timer.timeout.connect(self.start_build)
            self._json_timer = QTimer(self)
            self._json_timer.setSingleShot(True)
            self._json_timer.setInterval(500)
            self._json_timer.timeout.connect(self.apply_json_editor)
            self._export_timer = QTimer(self)
            self._export_timer.setSingleShot(True)
            self._export_timer.setInterval(200)
            self._export_timer.timeout.connect(self._start_export_build)
            self._setup_ui()
            self._load_initial_files()

        @staticmethod
        def _load_preferences() -> dict[str, Any]:
            """Load display preferences saved by this Windows user."""

            settings = QSettings("CCkit", "GrainBoundaryViewer")
            raw_pair_rules = str(settings.value("bond_pairs", "{}"))
            try:
                pair_rules = json.loads(raw_pair_rules) if raw_pair_rules.strip() else {}
            except json.JSONDecodeError:
                pair_rules = {}
            if not isinstance(pair_rules, dict):
                pair_rules = {}
            return {
                "atom_radius_scale": float(
                    settings.value("atom_radius_scale", 1.0)
                ),
                "show_bonds": str(settings.value("show_bonds", "false")).lower()
                in {"1", "true", "yes"},
                "bond_minimum": float(settings.value("bond_minimum", 0.1)),
                "bond_maximum": float(settings.value("bond_maximum", 2.5)),
                "bond_pairs": pair_rules,
                "show_plane": str(settings.value("show_plane", "true")).lower()
                in {"1", "true", "yes"},
                "plane_opacity": float(settings.value("plane_opacity", 0.25)),
                "plane_color": str(settings.value("plane_color", "#f59e0b")),
            }

        @staticmethod
        def _save_preferences(preferences: dict[str, Any]) -> None:
            """Persist display preferences in the Windows user profile."""

            settings = QSettings("CCkit", "GrainBoundaryViewer")
            for key, value in preferences.items():
                if key == "bond_pairs":
                    value = json.dumps(value, ensure_ascii=False)
                settings.setValue(key, value)
            settings.sync()

        @staticmethod
        def _dialog_directory(fallback: Path) -> Path:
            """Reuse the last selected folder across dialogs and app restarts."""

            settings = QSettings("CCkit", "GrainBoundaryViewer")
            saved = str(settings.value("files/last_directory", ""))
            if saved:
                directory = Path(saved).expanduser()
                if directory.is_dir():
                    return directory
            return fallback

        @staticmethod
        def _remember_file_directory(path: Path) -> None:
            """Persist the folder from an accepted file selection."""

            settings = QSettings("CCkit", "GrainBoundaryViewer")
            settings.setValue("files/last_directory", str(path.parent))
            settings.sync()

        def _setup_ui(self) -> None:
            root = QWidget(self)
            root_layout = QHBoxLayout(root)
            root_layout.setContentsMargins(8, 8, 8, 8)
            splitter = QSplitter(Qt.Horizontal, root)
            root_layout.addWidget(splitter)

            left = QWidget(splitter)
            left.setStyleSheet("background-color: #ffffff; color: #0f172a;")
            left_layout = QVBoxLayout(left)
            left_layout.setContentsMargins(0, 0, 0, 0)
            self.viewport = QtInteractor(left)
            self.viewport.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.viewport.set_background("#ffffff")
            self.viewport.render_window.SetAlphaBitPlanes(1)
            self.viewport.renderer.SetBackgroundAlpha(0.0)
            self.viewport.image_transparent_background = True
            left_layout.addWidget(self.viewport.interactor)
            frame_row = QHBoxLayout()
            frame_row.addWidget(QLabel("Frame"))
            self.frame_slider = QSlider(Qt.Orientation.Horizontal)
            self.frame_slider.setRange(0, 0)
            self.frame_slider.setEnabled(False)
            self.frame_slider.setTracking(False)
            self.frame_slider.setToolTip("Frame index starts at 0")
            self.frame_slider.valueChanged.connect(self.frame_changed)
            frame_row.addWidget(self.frame_slider, 1)
            self.frame_spin = QSpinBox()
            self.frame_spin.setRange(0, 0)
            self.frame_spin.setEnabled(False)
            self.frame_spin.setKeyboardTracking(False)
            self.frame_spin.valueChanged.connect(self.frame_changed)
            frame_row.addWidget(self.frame_spin)
            self.frame_count_label = QLabel("/ 0")
            frame_row.addWidget(self.frame_count_label)
            left_layout.addLayout(frame_row)
            self.render_status = QLabel("Select an input structure and gbin.json")
            self.render_status.setWordWrap(True)
            left_layout.addWidget(self.render_status)

            right_scroll = QScrollArea(splitter)
            right_scroll.setWidgetResizable(True)
            right_scroll.setMinimumWidth(430)
            panel = QWidget()
            panel_layout = QVBoxLayout(panel)
            panel_layout.setContentsMargins(10, 10, 10, 10)
            panel_layout.setSpacing(8)
            panel_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

            feature_group = QGroupBox("Visualization feature")
            feature_group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
            feature_layout = QFormLayout(feature_group)
            self.feature_combo = QComboBox()
            for feature in FEATURE_REGISTRY.all():
                self.feature_combo.addItem(feature.name, feature.key)
            self.feature_combo.currentIndexChanged.connect(self._feature_changed)
            feature_layout.addRow("Feature", self.feature_combo)
            self.feature_description = QLabel()
            self.feature_description.setWordWrap(True)
            self.feature_description.setStyleSheet("color: #475569;")
            feature_layout.addRow(self.feature_description)
            panel_layout.addWidget(feature_group, 0, Qt.AlignmentFlag.AlignTop)

            self.feature_stack = FeatureCardStack()
            self.omnigb_panel = QWidget()
            omnigb_layout = QVBoxLayout(self.omnigb_panel)
            omnigb_layout.setContentsMargins(0, 0, 0, 0)
            panel_layout.addWidget(self.feature_stack, 0, Qt.AlignmentFlag.AlignTop)

            files = QGroupBox("Input files")
            files_layout = QFormLayout(files)
            self.input_label = QLabel("not selected")
            self.input_label.setWordWrap(True)
            input_button = QPushButton("Choose structure…")
            input_button.clicked.connect(self.choose_input)
            files_layout.addRow(self.input_label, input_button)
            self.config_label = QLabel("not selected")
            self.config_label.setWordWrap(True)
            config_button = QPushButton("Choose gbin.json…")
            config_button.clicked.connect(self.choose_config)
            files_layout.addRow(self.config_label, config_button)
            omnigb_layout.addWidget(files)

            mode_group = QGroupBox("OmniGB mode")
            mode_layout = QFormLayout(mode_group)
            self.mode_combo = QComboBox()
            self.mode_combo.addItem("Build boundary", "build")
            self.mode_combo.addItem("Scan CSL candidates", "scan")
            self.mode_combo.currentIndexChanged.connect(self.mode_changed)
            mode_layout.addRow("Mode", self.mode_combo)
            omnigb_layout.addWidget(mode_group)

            common_group = QGroupBox("Common CSL parameters")
            common_layout = QVBoxLayout(common_group)
            axis_row = QHBoxLayout()
            axis_row.addWidget(QLabel("Rotation axis"))
            self.rotation_axis_editor = CrystalIndexEditor("direction")
            self.rotation_axis_editor.changed.connect(self.controls_changed)
            self.rotation_axis_editor.invalid.connect(self.interface_failed)
            axis_row.addWidget(self.rotation_axis_editor)
            common_layout.addLayout(axis_row)
            self.max_sigma_value = self._make_scalar_spin(
                common_layout,
                "Maximum Sigma",
                minimum=1.0,
                maximum=100000.0,
                decimals=0,
                step=1.0,
            )
            csl_row = QHBoxLayout()
            csl_row.addWidget(QLabel("CSL policy"))
            self.csl_mode_combo = QComboBox()
            self.csl_mode_combo.addItem("Near CSL", "near")
            self.csl_mode_combo.addItem("Exact CSL", "exact")
            self.csl_mode_combo.currentIndexChanged.connect(self.controls_changed)
            csl_row.addWidget(self.csl_mode_combo)
            common_layout.addLayout(csl_row)
            self.tolerance_values: dict[str, Any] = {}
            for name in ("du", "area", "deformation", "integer"):
                self.tolerance_values[name] = self._make_scalar_spin(
                    common_layout,
                    f"Tolerance {name}",
                    minimum=0.00001,
                    maximum=100.0,
                    decimals=5,
                    step=0.01,
                )
            omnigb_layout.addWidget(common_group)

            self.mode_stack = QStackedWidget()
            build_group = QGroupBox("Build parameters")
            build_layout = QVBoxLayout(build_group)
            self.build_angle_value, self.build_angle_step = self._make_value_row(
                build_layout, "Rotation angle (°)"
            )
            self.build_sigma_value = self._make_scalar_spin(
                build_layout,
                "Target Sigma",
                minimum=1.0,
                maximum=100000.0,
                decimals=0,
                step=1.0,
            )
            plane_row = QHBoxLayout()
            plane_row.addWidget(QLabel("Symmetric GB plane"))
            self.gb_plane_editor = CrystalIndexEditor("plane")
            self.gb_plane_editor.changed.connect(self.controls_changed)
            self.gb_plane_editor.invalid.connect(self.interface_failed)
            plane_row.addWidget(self.gb_plane_editor)
            build_layout.addLayout(plane_row)
            self.mode_stack.addWidget(build_group)

            scan_group = QGroupBox("Scan parameters")
            scan_layout = QVBoxLayout(scan_group)
            self.scan_start_value, self.scan_start_step = self._make_value_row(
                scan_layout, "Angle start (°)"
            )
            self.scan_stop_value, self.scan_stop_step = self._make_value_row(
                scan_layout, "Angle stop (°)"
            )
            self.scan_step_value, self.scan_step_step = self._make_value_row(
                scan_layout, "Angle increment (°)"
            )
            self.scan_step_value.setMinimum(0.00001)
            scan_hint = QLabel(
                "Only angles with a successful CSL match are retained in the CSV."
            )
            scan_hint.setWordWrap(True)
            scan_hint.setStyleSheet("color: #475569;")
            scan_layout.addWidget(scan_hint)
            self.scan_result_combo = QComboBox()
            self.scan_result_combo.setEnabled(False)
            scan_layout.addWidget(self.scan_result_combo)
            scan_actions = QHBoxLayout()
            self.run_scan_button = QPushButton("Run scan now")
            self.run_scan_button.clicked.connect(self.run_scan_now)
            self.run_scan_button.setEnabled(False)
            scan_actions.addWidget(self.run_scan_button)
            self.use_scan_result_button = QPushButton("Use selected in build")
            self.use_scan_result_button.clicked.connect(self.use_scan_result)
            self.use_scan_result_button.setEnabled(False)
            scan_actions.addWidget(self.use_scan_result_button)
            scan_layout.addLayout(scan_actions)
            self.scan_results_editor = QPlainTextEdit()
            self.scan_results_editor.setReadOnly(True)
            self.scan_results_editor.setPlaceholderText(
                "Successful scan rows will appear here."
            )
            self.scan_results_editor.setMinimumHeight(130)
            scan_layout.addWidget(self.scan_results_editor)
            self.export_scan_button = QPushButton("Export GB_scan_results.csv…")
            self.export_scan_button.clicked.connect(self.export_scan_results)
            self.export_scan_button.setEnabled(False)
            scan_layout.addWidget(self.export_scan_button)
            self.mode_stack.addWidget(scan_group)
            omnigb_layout.addWidget(self.mode_stack)

            summary = QGroupBox("Model summary")
            summary_layout = QFormLayout(summary)
            self.summary_mode = QLabel("—")
            self.summary_angle = QLabel("—")
            self.summary_sigma = QLabel("—")
            self.summary_axis = QLabel("—")
            self.summary_plane = QLabel("—")
            summary_layout.addRow("Mode", self.summary_mode)
            summary_layout.addRow("Angle", self.summary_angle)
            summary_layout.addRow("Sigma", self.summary_sigma)
            summary_layout.addRow("Axis", self.summary_axis)
            summary_layout.addRow("GB plane", self.summary_plane)
            omnigb_layout.addWidget(summary)

            self.control_group = QGroupBox("Live adjustment")
            control_layout = QVBoxLayout(self.control_group)
            self.dp1_value, self.dp1_step = self._make_value_row(
                control_layout, "dp1 (Å)"
            )
            self.dp2_value, self.dp2_step = self._make_value_row(
                control_layout, "dp2 (Å)"
            )
            self.cnid_grid_values = []
            for direction in ("u1", "u2"):
                row = QHBoxLayout()
                row.addWidget(QLabel(f"CNID {direction} subdivisions"), 1)
                spin = QSpinBox()
                spin.setRange(1, 1000)
                spin.setValue(4)
                spin.setKeyboardTracking(False)
                spin.valueChanged.connect(self.controls_changed)
                row.addWidget(spin)
                control_layout.addLayout(row)
                self.cnid_grid_values.append(spin)
            self.cnid_grid_label = QLabel()
            self.cnid_grid_label.setWordWrap(True)
            control_layout.addWidget(self.cnid_grid_label)
            self.normal_distance_value, self.normal_distance_step = self._make_value_row(
                control_layout, "GB spacing (Å)"
            )
            self.normal_distance_value.setMinimum(0.0)
            self.vacuum_value, self.vacuum_step = self._make_value_row(
                control_layout, "External vacuum (Å)"
            )
            self.vacuum_value.setMinimum(0.0)
            self.normal_ortho_check = QCheckBox(
                "Use an approximately normal GB cell vector"
            )
            self.normal_ortho_check.toggled.connect(self.controls_changed)
            control_layout.addWidget(self.normal_ortho_check)
            self.plane_ortho_check = QCheckBox(
                "Use approximately orthogonal in-plane vectors"
            )
            self.plane_ortho_check.toggled.connect(self.controls_changed)
            control_layout.addWidget(self.plane_ortho_check)
            self.mirror_check = QCheckBox(
                "Mirror the second grain (InterfaceMaster mirror=True)"
            )
            self.mirror_check.toggled.connect(self.controls_changed)
            control_layout.addWidget(self.mirror_check)
            omnigb_layout.addWidget(self.control_group)

            self.auto_check = QPushButton("Auto refresh: on")
            self.auto_check.setCheckable(True)
            self.auto_check.setChecked(True)
            self.auto_check.clicked.connect(self.toggle_auto_refresh)
            omnigb_layout.addWidget(self.auto_check)

            editor_group = QGroupBox("gbin.json")
            editor_layout = QVBoxLayout(editor_group)
            self.json_editor = QPlainTextEdit()
            self.json_editor.setLineWrapMode(QPlainTextEdit.NoWrap)
            self.json_editor.setMinimumHeight(320)
            self.json_editor.textChanged.connect(self.schedule_json_refresh)
            editor_layout.addWidget(self.json_editor)
            editor_buttons = QHBoxLayout()
            apply_button = QPushButton("Apply and refresh")
            apply_button.clicked.connect(self.apply_json_editor)
            save_button = QPushButton("Save gbin.json")
            save_button.clicked.connect(self.save_config)
            save_as_button = QPushButton("Save current JSON as…")
            save_as_button.clicked.connect(self.save_config_as)
            editor_buttons.addWidget(apply_button)
            editor_buttons.addWidget(save_button)
            editor_buttons.addWidget(save_as_button)
            editor_layout.addLayout(editor_buttons)
            omnigb_layout.addWidget(editor_group)

            self.feature_stack.addWidget(self.omnigb_panel)
            for feature in FEATURE_REGISTRY.all():
                if feature.key == "omnigb":
                    continue
                if feature.key == "wyckoff":
                    from core.wyckoff.wyckoff_panel import WyckoffPanel

                    self.wyckoff_panel = WyckoffPanel(
                        _load_model_module(feature.script_path), QSpinBox,
                        QDoubleSpinBox, QComboBox,
                        self._dialog_directory, self._remember_file_directory,
                    )
                    self.wyckoff_panel.changed.connect(self.wyckoff_changed)
                    self.wyckoff_panel.build_requested.connect(self.build_wyckoff_now)
                    self.wyckoff_panel.failed.connect(self.interface_failed)
                    self.feature_stack.addWidget(self.wyckoff_panel)
                    continue
                if feature.key == "interface":
                    from core.interface.interface_panel import InterfacePanel

                    self.interface_panel = InterfacePanel(
                        _load_model_module(feature.script_path), CrystalIndexEditor,
                        QSpinBox, QDoubleSpinBox, QComboBox, STRUCTURE_FILE_FILTER,
                        self._dialog_directory, self._remember_file_directory,
                    )
                    self.interface_panel.changed.connect(self.interface_changed)
                    self.interface_panel.build_requested.connect(self.build_interface_now)
                    self.interface_panel.failed.connect(self.interface_failed)
                    self.feature_stack.addWidget(self.interface_panel)
                    continue
                placeholder = QGroupBox(feature.card_title)
                placeholder_layout = QVBoxLayout(placeholder)
                placeholder_label = QLabel(
                    f"{feature.description}\n\nThis feature is not connected yet."
                )
                placeholder_label.setWordWrap(True)
                placeholder_layout.addWidget(placeholder_label)
                placeholder_layout.addStretch(1)
                self.feature_stack.addWidget(placeholder)

            self.error_label = QLabel()
            self.error_label.setWordWrap(True)
            self.error_label.setStyleSheet("color: #b42318;")
            panel_layout.addWidget(self.error_label)
            panel_layout.addStretch(1)
            right_scroll.setWidget(panel)
            self._feature_changed(0)
            splitter.setStretchFactor(0, 1)
            splitter.setStretchFactor(1, 0)
            self.setCentralWidget(root)

            toolbar = QToolBar("Personal preferences", self)
            toolbar.setMovable(False)
            load_structure_action = QAction("Load structure…", self)
            load_structure_action.triggered.connect(self.choose_input)
            toolbar.addAction(load_structure_action)
            load_trajectory_action = QAction("Load trajectory…", self)
            load_trajectory_action.triggered.connect(self.choose_trajectory)
            toolbar.addAction(load_trajectory_action)
            load_config_action = QAction("Load config…", self)
            load_config_action.triggered.connect(self.choose_config)
            toolbar.addAction(load_config_action)
            toolbar.addSeparator()
            preferences_action = QAction("Personal preferences", self)
            preferences_action.triggered.connect(self.open_preferences)
            toolbar.addAction(preferences_action)
            export_action = QAction("Export structure / trajectory…", self)
            export_action.setEnabled(False)
            export_action.triggered.connect(self.export_structure)
            toolbar.addAction(export_action)
            self.export_action = export_action
            self.addToolBar(Qt.TopToolBarArea, toolbar)

            file_menu = self.menuBar().addMenu("File")
            save_action = QAction("Save gbin.json", self)
            save_action.triggered.connect(self.save_config)
            file_menu.addAction(save_action)
            save_as_action = QAction("Save current JSON as…", self)
            save_as_action.triggered.connect(self.save_config_as)
            file_menu.addAction(save_as_action)
            export_menu_action = QAction("Export structure / trajectory…", self)
            export_menu_action.setEnabled(False)
            export_menu_action.triggered.connect(self.export_structure)
            file_menu.addAction(export_menu_action)
            self.export_menu_action = export_menu_action

        def open_preferences(self) -> None:
            """Open the top-level rendering preference dialog."""

            elements = (
                sorted({str(site.specie) for site in self.last_structure})
                if self.last_structure is not None
                else []
            )
            dialog = PreferenceDialog(self.preferences, elements, self)
            if dialog.exec() != QDialog.DialogCode.Accepted:
                return
            self.preferences = dialog.values()
            self._save_preferences(self.preferences)
            if self.last_structure is not None:
                self._render_structure(self.last_structure)

        def _make_value_row(self, layout: Any, label: str) -> tuple[Any, Any]:
            row = QHBoxLayout()
            value_spin = QDoubleSpinBox()
            value_spin.setRange(-10000.0, 10000.0)
            value_spin.setDecimals(5)
            value_spin.setSingleStep(0.1)
            value_spin.valueChanged.connect(self.controls_changed)
            step_spin = QDoubleSpinBox()
            step_spin.setRange(0.00001, 10000.0)
            step_spin.setDecimals(5)
            step_spin.setValue(0.1)
            step_spin.valueChanged.connect(value_spin.setSingleStep)
            row.addWidget(value_spin, 3)
            row.addWidget(QLabel("step"))
            row.addWidget(step_spin, 2)
            layout.addLayout(self._labelled_row(label, row))
            return value_spin, step_spin

        def _make_scalar_spin(
            self,
            layout: Any,
            label: str,
            *,
            minimum: float,
            maximum: float,
            decimals: int,
            step: float,
        ) -> Any:
            """Create a labeled scalar editor used by both OmniGB modes."""

            spin = QDoubleSpinBox()
            spin.setRange(minimum, maximum)
            spin.setDecimals(decimals)
            spin.setSingleStep(step)
            spin.valueChanged.connect(self.controls_changed)
            layout.addLayout(self._labelled_row(label, self._spin_layout(spin)))
            return spin

        @staticmethod
        def _spin_layout(spin: Any) -> Any:
            row = QHBoxLayout()
            row.addWidget(spin)
            return row

        def _make_vector_row(
            self, layout: Any, label: str, *, integer: bool
        ) -> list[Any]:
            """Create three Cartesian component editors for an axis or plane."""

            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            values = []
            for component in "xyz":
                row.addWidget(QLabel(component))
                spin = QDoubleSpinBox()
                spin.setRange(-10000.0, 10000.0)
                spin.setDecimals(0 if integer else 5)
                spin.setSingleStep(1.0 if integer else 0.1)
                spin.valueChanged.connect(self.controls_changed)
                row.addWidget(spin)
                values.append(spin)
            layout.addLayout(row)
            return values

        @staticmethod
        def _labelled_row(label: str, row: Any) -> Any:
            container = QHBoxLayout()
            container.addWidget(QLabel(label))
            container.addLayout(row)
            return container

        def _active_input_path(self) -> Path | None:
            if self.active_feature.key == "wyckoff":
                return self.wyckoff_panel.input_path
            return (
                self.interface_panel.input_path if self.active_feature.key == "interface"
                else self.input_path
            )

        def _active_auto_refresh(self) -> bool:
            if self.active_feature.key == "wyckoff":
                return False
            return (
                self.interface_panel.auto_check.isChecked()
                if self.active_feature.key == "interface" else self.auto_refresh
            )

        def _active_config_text(self) -> str:
            if self.active_feature.key == "wyckoff":
                return json.dumps(self.wyckoff_panel.config())
            return (
                self.interface_panel.editor.toPlainText()
                if self.active_feature.key == "interface" else self.json_editor.toPlainText()
            )

        def interface_failed(self, error: str) -> None:
            if hasattr(self, "error_label"):
                self.error_label.setText(error)

        def interface_changed(self) -> None:
            if self.active_feature.key == "interface":
                self.queue_build()

        def build_interface_now(self) -> None:
            if self.active_feature.key != "interface":
                return
            self.queue_build()
            if getattr(self, "_pending_config", None) is not None:
                self._build_timer.start()

        def wyckoff_changed(self) -> None:
            """Invalidate an older request while keeping generation manual."""

            self.feature_stack.updateGeometry()
            if self.active_feature.key == "wyckoff":
                self._build_timer.stop()
                self._pending_config = None
                self.build_serial += 1
                if self.worker is not None and self.worker.isRunning() and not self._export_build_active:
                    self.worker.requestInterruption()

        def build_wyckoff_now(self) -> None:
            if self.active_feature.key != "wyckoff":
                return
            self.queue_build()
            if getattr(self, "_pending_config", None) is not None:
                self._build_timer.start()

        def _feature_changed(self, index: int) -> None:
            """Switch independently owned cards and invalidate pending old builds."""

            self._build_timer.stop()
            self._json_timer.stop()
            self._pending_config = None
            self.build_serial += 1
            if self.worker is not None and self.worker.isRunning() and not self._export_build_active:
                self.worker.requestInterruption()
            feature = FEATURE_REGISTRY.all()[index]
            self.active_feature = feature
            self.feature_stack.setCurrentIndex(index)
            self.feature_description.setText(feature.description)
            if not feature.available:
                self.render_status.setText(f"{feature.name} is reserved for a future CCkit script.")
            elif feature.key == "wyckoff":
                self.error_label.clear()
                self.render_status.setText("Fill in Wyckoff parameters, then click Generate crystal")
            else:
                self.queue_build()

        def select_feature(self, key: str) -> None:
            """Select a feature by stable registry key for future integrations."""

            for index, feature in enumerate(FEATURE_REGISTRY.all()):
                if feature.key == key:
                    self.feature_combo.setCurrentIndex(index)
                    return
            raise KeyError(f"unknown viewer feature: {key}")

        def _load_initial_files(self) -> None:
            if self.input_path is not None and self.input_path.is_file():
                self.input_label.setText(str(self.input_path))
            if self.config_path is not None and self.config_path.is_file():
                self.config_label.setText(str(self.config_path))
                try:
                    config = _read_config(self.config_path)
                except Exception as exc:
                    self.error_label.setText(str(exc))
                else:
                    self._set_editor_config(config)
                    self.queue_build()

        def choose_input(self) -> None:
            if self.active_feature.key == "wyckoff":
                self.render_status.setText("Wyckoff model uses the space group, lattice and sites in its card")
                return
            if self.active_feature.key == "interface":
                self.interface_panel.choose_structure("substrate")
                return
            directory = self._dialog_directory(
                self.input_path.parent if self.input_path else REPOSITORY_ROOT
            )
            selected, _ = QFileDialog.getOpenFileName(
                self,
                "Choose crystal structure",
                str(directory),
                STRUCTURE_FILE_FILTER,
            )
            if selected:
                self.input_path = Path(selected).resolve()
                self._remember_file_directory(self.input_path)
                self.input_label.setText(str(self.input_path))
                self.queue_build()

        def choose_trajectory(self) -> None:
            """Read a saved lattice trajectory without rebuilding it as bulk."""

            if self._export_request is not None:
                self.error_label.setText("Wait for the current export to finish.")
                return
            directory = self._dialog_directory(
                self.input_path.parent if self.input_path else REPOSITORY_ROOT
            )
            selected, _ = QFileDialog.getOpenFileName(
                self, "Open lattice trajectory", str(directory),
                "Trajectory files (*.xyz *.extxyz *.traj);;All files (*)",
            )
            if not selected:
                return
            path = Path(selected).expanduser().resolve()
            try:
                trajectory = _read_trajectory(path)
            except Exception as exc:
                self.error_label.setText(f"Could not load trajectory: {exc}")
                return
            self._build_timer.stop()
            self._json_timer.stop()
            self.build_serial += 1
            if self.worker is not None and self.worker.isRunning():
                self.worker.requestInterruption()
            self.trajectory_source = path
            self.rendered_config = None
            self.rendered_input_path = None
            self.rendered_feature_key = None
            self.error_label.clear()
            self._remember_file_directory(path)
            self._set_trajectory(trajectory)
            self._set_export_enabled(True)

        def choose_config(self) -> None:
            if self.active_feature.key == "wyckoff":
                self.wyckoff_panel.load_config()
                return
            if self.active_feature.key == "interface":
                self.interface_panel.load_config()
                return
            directory = self._dialog_directory(
                self.config_path.parent if self.config_path else REPOSITORY_ROOT
            )
            selected, _ = QFileDialog.getOpenFileName(
                self,
                "Choose gbin.json",
                str(directory),
                "JSON files (*.json);;All files (*)",
            )
            if selected:
                self.config_path = Path(selected).resolve()
                self._remember_file_directory(self.config_path)
                self.config_label.setText(str(self.config_path))
                try:
                    self._set_editor_config(_read_config(self.config_path))
                except Exception as exc:
                    self.error_label.setText(str(exc))
                else:
                    self.queue_build()

        def _set_editor_config(self, config: dict[str, Any]) -> None:
            config = _viewer_config(config)
            self.updating_controls = True
            try:
                self.json_editor.setPlainText(
                    json.dumps(config, indent=2, ensure_ascii=False)
                )
                self._sync_controls(config)
                self._update_summary(config)
            finally:
                self.updating_controls = False

        def _sync_controls(self, config: dict[str, Any]) -> None:
            try:
                model_values = _omnigb_parameters(config)
                dp1, dp2 = _termination_values(config)
                grid = _load_model_module().cnid_grid(config)
                normal_distance = float(
                    config.get("normal_distance", config.get("dx", 0.0))
                )
                vacuum = float(config.get("vacuum", config.get("vx", 0.0)))
                if normal_distance < 0 or vacuum < 0:
                    raise ValueError("normal_distance and vacuum must be nonnegative")
                normal_ortho = config.get("normal_ortho", True)
                plane_ortho = config.get("plane_ortho", False)
                if not isinstance(normal_ortho, bool):
                    raise ValueError("normal_ortho must be a JSON boolean")
                if not isinstance(plane_ortho, bool):
                    raise ValueError("plane_ortho must be a JSON boolean")
                mirror = config.get("mirror", False)
                if not isinstance(mirror, bool):
                    raise ValueError("mirror must be a JSON boolean")
            except (TypeError, ValueError) as exc:
                self.error_label.setText(str(exc))
                return
            for spin, value in (
                (self.dp1_value, dp1),
                (self.dp2_value, dp2),
                *zip(self.cnid_grid_values, grid),
                (self.normal_distance_value, normal_distance),
                (self.vacuum_value, vacuum),
                (self.max_sigma_value, model_values["max_sigma"]),
                (self.build_angle_value, model_values["rotation_angle_deg"]),
                (self.build_sigma_value, model_values["sigma"]),
                (self.scan_start_value, model_values["angle_start_deg"]),
                (self.scan_stop_value, model_values["angle_stop_deg"]),
                (self.scan_step_value, model_values["angle_step_deg"]),
            ):
                blocker = QSignalBlocker(spin)
                spin.setValue(value)
                del blocker
            self.rotation_axis_editor.set_indices(
                model_values["rotation_axis"], model_values["rotation_axis_notation"]
            )
            self.gb_plane_editor.set_indices(
                model_values["gb_plane"], model_values["gb_plane_notation"]
            )
            mode_blocker = QSignalBlocker(self.mode_combo)
            mode_index = self.mode_combo.findData(model_values["mode"])
            if mode_index >= 0:
                self.mode_combo.setCurrentIndex(mode_index)
            del mode_blocker
            csl_blocker = QSignalBlocker(self.csl_mode_combo)
            csl_index = self.csl_mode_combo.findData(model_values["csl_mode"])
            if csl_index >= 0:
                self.csl_mode_combo.setCurrentIndex(csl_index)
            del csl_blocker
            for name, value in model_values["tolerances"].items():
                blocker = QSignalBlocker(self.tolerance_values[name])
                self.tolerance_values[name].setValue(value)
                del blocker
            blocker = QSignalBlocker(self.mirror_check)
            self.mirror_check.setChecked(mirror)
            del blocker
            blocker = QSignalBlocker(self.normal_ortho_check)
            self.normal_ortho_check.setChecked(normal_ortho)
            del blocker
            blocker = QSignalBlocker(self.plane_ortho_check)
            self.plane_ortho_check.setChecked(plane_ortho)
            del blocker
            self.cnid_grid_label.setText(
                f"{grid[0] * grid[1]} frames · indices 0–{grid[0] * grid[1] - 1}"
            )
            self._set_mode_card(model_values["mode"])

        def _set_mode_card(self, mode: str) -> None:
            """Show the selected mode card without changing its saved values."""

            is_scan = mode == "scan"
            self.mode_stack.setCurrentIndex(1 if is_scan else 0)
            self.run_scan_button.setEnabled(is_scan)
            has_results = bool(self.scan_rows)
            self.scan_result_combo.setEnabled(is_scan and has_results)
            self.use_scan_result_button.setEnabled(is_scan and has_results)

        def mode_changed(self, _index: int) -> None:
            """Change only ``mode`` so all shared JSON parameters survive."""

            if self.updating_controls:
                return
            mode = str(self.mode_combo.currentData() or "build")
            try:
                config = json.loads(self.json_editor.toPlainText())
                if not isinstance(config, dict):
                    raise ValueError("gbin.json must contain a JSON object")
                config["mode"] = mode
                self.updating_controls = True
                try:
                    self.json_editor.setPlainText(
                        json.dumps(config, indent=2, ensure_ascii=False)
                    )
                finally:
                    self.updating_controls = False
                self._set_mode_card(mode)
                self._update_summary(config)
                self.error_label.clear()
            except (json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"JSON/config error: {exc}")
                return
            self.build_serial += 1
            if self.auto_refresh:
                self.queue_build(config)

        def _update_summary(self, config: dict[str, Any]) -> None:
            mode = str(config.get("mode", "build"))
            self.summary_mode.setText(mode)
            if mode == "scan":
                start = config.get("angle_start_deg", 0.0)
                stop = config.get("angle_stop_deg", 180.0)
                self.summary_angle.setText(f"{start}° → {stop}°")
                self.summary_sigma.setText(
                    f"≤ {config.get('max_sigma', config.get('maxsigma', '—'))}"
                )
            else:
                self.summary_angle.setText(
                    f"{config.get('rotation_angle_deg', config.get('angle_deg', '—'))}°"
                )
                self.summary_sigma.setText(str(config.get("sigma", "—")))
            self.summary_axis.setText(
                str(config.get("rotation_axis", config.get("axis", "—")))
            )
            self.summary_plane.setText(
                str(config.get("gb_plane", config.get("plane", "—")))
            )

        def schedule_json_refresh(self) -> None:
            if not self.updating_controls:
                self._json_timer.start()

        def apply_json_editor(self) -> None:
            try:
                config = json.loads(self.json_editor.toPlainText())
                if not isinstance(config, dict):
                    raise ValueError("gbin.json must contain a JSON object")
                if config.get("mode", "build") not in {"build", "scan"}:
                    raise ValueError("mode must be either 'build' or 'scan'")
                self.error_label.clear()
                config = _viewer_config(config)
                self._set_editor_config(config)
            except (json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"JSON/config error: {exc}")
                return
            if self.auto_refresh:
                self.queue_build(config)

        def _update_config_from_controls(self, config: dict[str, Any]) -> dict[str, Any]:
            """Copy every visible OmniGB control into one JSON object."""

            shifts = config.setdefault("termination_shift", {})
            if not isinstance(shifts, dict):
                raise ValueError("termination_shift must be a JSON object")
            if "dp1" in config:
                config["dp1"] = self.dp1_value.value()
            else:
                shifts["dp1"] = self.dp1_value.value()
            if "dp2" in config:
                config["dp2"] = self.dp2_value.value()
            else:
                shifts["dp2"] = self.dp2_value.value()
            config.pop("translation", None)
            config.pop("dydz", None)
            config["cnid_grid"] = [spin.value() for spin in self.cnid_grid_values]
            count = config["cnid_grid"][0] * config["cnid_grid"][1]
            self.cnid_grid_label.setText(f"{count} frames · indices 0–{count - 1}")
            config["normal_distance"] = self.normal_distance_value.value()
            config["vacuum"] = self.vacuum_value.value()
            config["normal_ortho"] = self.normal_ortho_check.isChecked()
            config["plane_ortho"] = self.plane_ortho_check.isChecked()
            config["mirror"] = self.mirror_check.isChecked()
            config["rotation_axis"] = self.rotation_axis_editor.indices()
            config["rotation_axis_notation"] = self.rotation_axis_editor.notation()
            config["max_sigma"] = int(round(self.max_sigma_value.value()))
            config["csl_mode"] = str(self.csl_mode_combo.currentData() or "near")
            if "exact" in config:
                config["exact"] = config["csl_mode"] == "exact"
            config["tolerances"] = {
                name: spin.value() for name, spin in self.tolerance_values.items()
            }
            config["rotation_angle_deg"] = self.build_angle_value.value()
            config["sigma"] = int(round(self.build_sigma_value.value()))
            config["gb_plane"] = self.gb_plane_editor.indices()
            config["gb_plane_notation"] = self.gb_plane_editor.notation()
            config["angle_start_deg"] = self.scan_start_value.value()
            config["angle_stop_deg"] = self.scan_stop_value.value()
            config["angle_step_deg"] = self.scan_step_value.value()
            return config

        def controls_changed(self, *_args: Any) -> None:
            if self.updating_controls or not hasattr(self, "json_editor"):
                return
            try:
                config = json.loads(self.json_editor.toPlainText())
                if not isinstance(config, dict):
                    raise ValueError("gbin.json must contain a JSON object")
                if config.get("mode", "build") not in {"build", "scan"}:
                    raise ValueError("mode must be either 'build' or 'scan'")
                config = self._update_config_from_controls(config)
                self.updating_controls = True
                try:
                    self.json_editor.setPlainText(
                        json.dumps(config, indent=2, ensure_ascii=False)
                    )
                finally:
                    self.updating_controls = False
                self._set_mode_card(str(config.get("mode", "build")))
                self._update_summary(config)
                self.error_label.clear()
            except (json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"JSON/config error: {exc}")
                return
            if self.auto_refresh:
                self.queue_build(config)

        def toggle_auto_refresh(self, checked: bool) -> None:
            self.auto_refresh = checked
            self.auto_check.setText(
                f"Auto refresh: {'on' if checked else 'off'}"
            )
            if checked:
                self.apply_json_editor()

        def queue_build(self, config: dict[str, Any] | None = None) -> None:
            if not self.active_feature.available:
                return
            source = self._active_input_path()
            if source is None or (self.active_feature.key != "wyckoff" and not source.is_file()):
                self._pending_config = None
                if self.active_feature.key == "interface":
                    self.render_status.setText("Select substrate and film structures")
                return
            try:
                if config is None:
                    text = self._active_config_text()
                    if not text.strip():
                        self._pending_config = None
                        return
                    config = json.loads(text)
                if not isinstance(config, dict):
                    raise ValueError("configuration must be a JSON object")
                config = (
                    self.wyckoff_panel.backend.parameters(config)
                    if self.active_feature.key == "wyckoff" else _viewer_config(config)
                )
                if self.active_feature.key == "interface":
                    config = self.interface_panel.backend.parameters(config)
                    if not config.get("film_structure"):
                        self._pending_config = None
                        self.render_status.setText("Select a film structure")
                        return
            except (OSError, RuntimeError, TypeError, ValueError) as exc:
                self._pending_config = None
                self.error_label.setText(f"JSON/config error: {exc}")
                return
            self._pending_config = copy.deepcopy(config)
            if self.worker is not None and self.worker.isRunning() and not self._export_build_active:
                self.worker.requestInterruption()
            self.build_serial += 1
            if self._active_auto_refresh():
                self._build_timer.start()

        def start_build(self) -> None:
            if self._export_request is not None:
                self._start_export_build()
                return
            source = self._active_input_path()
            if source is None or (self.active_feature.key != "wyckoff" and not source.is_file()):
                return
            config = getattr(self, "_pending_config", None)
            if not isinstance(config, dict):
                return
            if (
                self.worker is not None
                and self.worker.isRunning()
            ) or (
                self.scan_worker is not None
                and self.scan_worker.isRunning()
            ):
                # Do not start two InterfaceMaster operations concurrently.
                # The timer will pick up the newest pending JSON after the
                # current build finishes.
                self._build_timer.start(200)
                return
            if config.get("mode", "build") == "scan":
                self.start_scan(config)
                return
            self.build_serial += 1
            serial = self.build_serial
            self.render_status.setText("Building preview…")
            self.worker = BuildWorker(source, config, self.active_feature)
            worker = self.worker
            self.worker.completed.connect(
                lambda structure, token=serial, settings=worker.config,
                source=worker.input_path, key=worker.feature.key: self.preview_ready(
                    token, structure, settings, source, key
                )
            )
            self.worker.failed.connect(
                lambda error, token=serial: self.preview_failed(token, error)
            )
            self.worker.progress.connect(
                lambda completed, total, token=serial: self.build_progress(token, completed, total)
            )
            self.worker.start()

        def build_progress(self, serial: int, completed: int, total: int) -> None:
            """Show frame generation progress only for the active preview."""

            if serial == self.build_serial:
                message = (
                    f"Expanding Wyckoff positions: {completed}/{total}…"
                    if self.active_feature.key == "wyckoff"
                    else "Matching contact planes…"
                    if completed == 0 and self.active_feature.key == "interface"
                    else f"Sampling CNID: {completed}/{total} frames…"
                )
                self.render_status.setText(message)

        def start_scan(self, config: dict[str, Any] | None = None) -> None:
            """Start a scan for the current JSON configuration."""

            if self.input_path is None or not self.input_path.is_file():
                self.error_label.setText("Choose an input structure before scanning.")
                return
            if (
                self.worker is not None
                and self.worker.isRunning()
            ) or (
                self.scan_worker is not None
                and self.scan_worker.isRunning()
            ):
                self._pending_config = config or getattr(
                    self, "_pending_config", None
                )
                self.build_serial += 1
                self._build_timer.start(200)
                return
            config = config or getattr(self, "_pending_config", None)
            if not isinstance(config, dict):
                return
            self.build_serial += 1
            serial = self.build_serial
            self.render_status.setText("Scanning CSL candidates…")
            self.scan_rows = []
            self.scan_csv_text = ""
            self.scan_result_combo.clear()
            self.scan_results_editor.clear()
            self.export_scan_button.setEnabled(False)
            self.scan_worker = ScanWorker(
                self.input_path, config, self.active_feature
            )
            self.scan_worker.completed.connect(
                lambda rows, csv_text, token=serial: self.scan_ready(
                    token, rows, csv_text
                )
            )
            self.scan_worker.failed.connect(
                lambda error, token=serial: self.scan_failed(token, error)
            )
            self.scan_worker.start()

        def run_scan_now(self) -> None:
            """Run scan mode immediately, even when auto refresh is disabled."""

            try:
                config = json.loads(self.json_editor.toPlainText())
                if not isinstance(config, dict):
                    raise ValueError("gbin.json must contain a JSON object")
                config["mode"] = "scan"
                self._pending_config = config
            except (json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"JSON/config error: {exc}")
                return
            self.start_scan(config)

        def scan_failed(self, serial: int, error: str) -> None:
            if serial != self.build_serial:
                return
            self.render_status.setText("CSL scan failed")
            self.error_label.setText(error)

        def scan_ready(
            self, serial: int, rows: list[dict[str, str]], csv_text: str
        ) -> None:
            if serial != self.build_serial:
                return
            try:
                current_config = json.loads(self.json_editor.toPlainText())
            except (json.JSONDecodeError, TypeError):
                return
            if not isinstance(current_config, dict) or current_config.get(
                "mode", "build"
            ) != "scan":
                return
            self.scan_rows = rows
            self.scan_csv_text = csv_text
            self.scan_results_editor.setPlainText(_format_scan_table(rows))
            self.scan_result_combo.clear()
            for row in rows:
                self.scan_result_combo.addItem(
                    "Sigma ({sigma_1}, {sigma_2}) · {angle}° · {mode}".format(
                        sigma_1=row.get("sigma_1", "?"),
                        sigma_2=row.get("sigma_2", "?"),
                        angle=row.get("angle_deg", "?"),
                        mode=row.get("csl_mode", "?"),
                    )
                )
            self._set_mode_card("scan")
            self.export_scan_button.setEnabled(bool(rows))
            self.render_status.setText(
                f"CSL scan complete · {len(rows)} successful candidates"
            )
            self.error_label.clear()

        def use_scan_result(self) -> None:
            """Copy one successful scan row into the build controls."""

            index = self.scan_result_combo.currentIndex()
            if index < 0 or index >= len(self.scan_rows):
                return
            row = self.scan_rows[index]
            try:
                angle = float(row["angle_deg"])
                sigma = int(row["sigma_1"])
                if sigma != int(row["sigma_2"]):
                    raise ValueError("the selected scan row has unequal Sigma values")
                config = json.loads(self.json_editor.toPlainText())
                if not isinstance(config, dict):
                    raise ValueError("gbin.json must contain a JSON object")
                config["mode"] = "build"
                config["rotation_angle_deg"] = angle
                config["sigma"] = sigma
                self._set_editor_config(config)
            except (KeyError, json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"Could not select scan result: {exc}")
                return
            if self.auto_refresh:
                self.queue_build(config)

        def export_scan_results(self) -> None:
            """Write the most recent successful scan rows to a CSV file."""

            if not self.scan_csv_text:
                self.error_label.setText("Run a scan before exporting its CSV.")
                return
            initial_directory = (
                self.input_path.parent if self.input_path else REPOSITORY_ROOT
            )
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Export GB scan results",
                str(self._dialog_directory(initial_directory) / "GB_scan_results.csv"),
                "CSV files (*.csv);;All files (*)",
            )
            if not selected:
                return
            target = Path(selected).expanduser().resolve()
            self._remember_file_directory(target)
            try:
                target.write_text(self.scan_csv_text, encoding="utf-8")
            except OSError as exc:
                self.error_label.setText(f"Could not export scan CSV: {exc}")
                return
            self.error_label.setText(f"Exported scan CSV: {target}")

        def preview_failed(self, serial: int, error: str) -> None:
            if serial != self.build_serial:
                return
            try:
                current_config = json.loads(self._active_config_text())
            except (RuntimeError, ValueError, TypeError):
                return
            if not isinstance(current_config, dict) or current_config.get(
                "mode", "build"
            ) != "build":
                return
            self.render_status.setText("Preview build failed")
            self.error_label.setText(error)

        def _request_matches_editor(
            self, config: dict[str, Any], input_path: Path, feature_key: str
        ) -> bool:
            """Compare a build with the currently selected card and its two inputs."""

            if self.active_feature.key != feature_key:
                return False
            try:
                current = json.loads(self._active_config_text())
                if not isinstance(current, dict):
                    return False
                current = (
                    self.wyckoff_panel.backend.parameters(current)
                    if feature_key == "wyckoff" else _viewer_config(current)
                )
                if feature_key == "interface":
                    current = self.interface_panel.backend.parameters(current)
            except (OSError, RuntimeError, TypeError, ValueError):
                return False
            return current == config and self._active_input_path() == input_path

        def preview_ready(
            self,
            serial: int,
            structure: Any,
            config: dict[str, Any],
            input_path: Path,
            feature_key: str,
        ) -> None:
            if serial != self.build_serial:
                return
            if (
                config.get("mode", "build") != "build"
                or not self._request_matches_editor(config, input_path, feature_key)
            ):
                return
            self.error_label.clear()
            self.trajectory_source = None
            self.rendered_config = copy.deepcopy(config)
            self.rendered_input_path = input_path
            self.rendered_feature_key = feature_key
            self._set_export_enabled(self._export_request is None)
            self._set_trajectory(structure)

        def _set_trajectory(self, trajectory: StructureTrajectory) -> None:
            """Install a sequence atomically and reset its frame index to zero."""

            if not trajectory.frames:
                raise ValueError("a trajectory must contain at least one frame")
            self.trajectory = trajectory
            if self.rendered_feature_key == "wyckoff" and trajectory.details:
                details = trajectory.details
                self.wyckoff_panel.summary.setText(
                    f"Space group {details['space_group']} · {details['symbol']} · "
                    f"{details['atom_count']} atoms"
                )
            if self.active_feature.key == "interface" and trajectory.details:
                details = trajectory.details
                self.interface_panel.summary.setText(
                    f"Twist {details['twist_angle_deg']:.4f}° · "
                    f"surface indices {details['substrate_area_index']}/{details['film_area_index']}\n"
                    f"Area {details['interface_area_A2']:.3f} Å² · "
                    f"film metric strain {details['film_metric_strain_max']:.5f} · "
                    f"{details['atom_count']} atoms"
                )
            maximum = len(trajectory.frames) - 1
            for control in (self.frame_slider, self.frame_spin):
                blocker = QSignalBlocker(control)
                control.setRange(0, maximum)
                control.setValue(0)
                control.setEnabled(maximum > 0)
                del blocker
            self.frame_count_label.setText(f"/ {maximum}")
            self.frame_count_label.setFixedWidth(
                self.frame_count_label.fontMetrics().horizontalAdvance(f"/ {maximum}") + 8
            )
            self.frame_changed(0)

        def frame_changed(self, index: int) -> None:
            """Synchronize slider and editor, rendering without rebuilding CSL."""

            if self.trajectory is None or not 0 <= index < len(self.trajectory.frames):
                return
            self.current_frame_index = index
            for control in (self.frame_slider, self.frame_spin):
                blocker = QSignalBlocker(control)
                control.setValue(index)
                del blocker
            structure = self.trajectory.frames[index]
            self.last_structure = structure
            self._render_structure(structure)
            metadata = self.trajectory.metadata[index]
            rbt = ""
            if "rbt_alpha" in metadata and "rbt_beta" in metadata:
                try:
                    alpha = float(metadata["rbt_alpha"])
                    beta = float(metadata["rbt_beta"])
                except (TypeError, ValueError):
                    pass
                else:
                    rbt = f" · RBT α={alpha:.4f}, β={beta:.4f}"
            self.render_status.setText(
                f"Frame {index}/{len(self.trajectory.frames) - 1}"
                f" · {len(self.trajectory.frames)} frames · {len(structure)} atoms{rbt}"
                f" · cell {structure.lattice.a:.2f} × {structure.lattice.b:.2f}"
                f" × {structure.lattice.c:.2f} Å"
            )

        def _render_structure(self, structure: Any) -> None:
            saved_camera: dict[str, Any] | None = None
            if self._camera_initialized:
                camera = self.viewport.camera
                saved_camera = {
                    "position": tuple(camera.position),
                    "focal_point": tuple(camera.focal_point),
                    "up": tuple(camera.up),
                    "parallel_scale": float(camera.parallel_scale),
                    "parallel_projection": bool(camera.parallel_projection),
                    "view_angle": float(camera.view_angle),
                }
            self.viewport.clear()
            points = np.asarray(structure.cart_coords, dtype=float)
            elements = [str(site.specie) for site in structure]
            palette = (
                "#4cc9f0",
                "#f72585",
                "#fca311",
                "#80ed99",
                "#b8c0ff",
                "#ffd166",
            )
            for index, element in enumerate(dict.fromkeys(elements)):
                selected = np.asarray([item == element for item in elements])
                cloud = pv.PolyData(points[selected])
                radius = max(
                    min(structure.lattice.abc)
                    * 0.025
                    * self.preferences["atom_radius_scale"],
                    0.12,
                )
                glyph = cloud.glyph(
                    geom=pv.Sphere(radius=radius),
                    scale=False,
                    orient=False,
                )
                self.viewport.add_mesh(
                    glyph,
                    color=palette[index % len(palette)],
                    smooth_shading=True,
                    name=f"atoms_{element}",
                    reset_camera=False,
                )
            if self.preferences["show_bonds"]:
                bond_mesh = _bond_lines(
                    np,
                    pv,
                    structure,
                    self.preferences["bond_minimum"],
                    self.preferences["bond_maximum"],
                    self.preferences.get("bond_pairs"),
                )
                if bond_mesh.n_points:
                    self.viewport.add_mesh(
                        bond_mesh,
                        color="#94a3b8",
                        line_width=2.0,
                        name="bonds",
                        reset_camera=False,
                    )
            if self.preferences["show_plane"]:
                plane_color = self.preferences["plane_color"]
                fractions = (0.0, 0.5)
                if self.trajectory is not None:
                    metadata = self.trajectory.metadata[self.current_frame_index]
                    fractions = metadata.get("interface_plane_fractions", fractions)
                planes = _grain_boundary_planes(np, pv, structure, fractions)
                for index, plane in enumerate(planes):
                    self.viewport.add_mesh(
                        plane,
                        color=plane_color,
                        opacity=self.preferences["plane_opacity"],
                        show_edges=True,
                        edge_color=plane_color,
                        name=f"gb_plane_{index}",
                        reset_camera=False,
                    )
                if fractions:
                    self.viewport.add_text(
                        "Interface planes at fractional x: " + ", ".join(f"{x:.4f}" for x in fractions),
                        position="upper_left",
                        font_size=11,
                        color=plane_color,
                        name="gb_plane_label",
                    )
            self.viewport.add_mesh(
                _lattice_lines(np, pv, structure.lattice),
                color="#475569",
                line_width=2.0,
                name="unit_cell",
                reset_camera=False,
            )
            axes = self.viewport.add_axes(
                xlabel="x",
                ylabel="y",
                zlabel="z",
                labels_off=False,
                color="#0f172a",
                x_color="#dc2626",
                y_color="#15803d",
                z_color="#2563eb",
                line_width=4,
                label_size=(0.45, 0.18),
                viewport=(0.0, 0.0, 0.25, 0.25),
            )
            for caption in (
                axes.GetXAxisCaptionActor2D(),
                axes.GetYAxisCaptionActor2D(),
                axes.GetZAxisCaptionActor2D(),
            ):
                text_property = caption.GetCaptionTextProperty()
                text_property.SetFontSize(24)
                text_property.SetBold(True)
                text_property.SetShadow(False)
            self.viewport.show_grid(
                color="#0f172a",
                location="outer",
                ticks="both",
                xtitle="x (Å)",
                ytitle="y (Å)",
                ztitle="z (Å)",
                font_size=18,
                bold=True,
                use_3d_text=False,
            )
            if saved_camera is None:
                self.viewport.reset_camera()
                self._camera_initialized = True
            else:
                camera = self.viewport.camera
                camera.position = saved_camera["position"]
                camera.focal_point = saved_camera["focal_point"]
                camera.up = saved_camera["up"]
                camera.parallel_scale = saved_camera["parallel_scale"]
                camera.parallel_projection = saved_camera["parallel_projection"]
                camera.view_angle = saved_camera["view_angle"]
                self.viewport.renderer.reset_camera_clipping_range()
            self.viewport.render()

        def export_structure(self) -> None:
            """Capture parameters and selected frame before opening export dialogs."""

            if self._export_request is not None:
                self.error_label.setText("The current export is still being built.")
                return
            self.frame_spin.interpretText()
            frame_index = self.current_frame_index
            cached = self.trajectory
            source = self.trajectory_source
            feature = self.active_feature
            try:
                if source is not None and cached is not None:
                    config = {}
                    input_path = source
                else:
                    config = copy.deepcopy(self._read_editor_config())
                    if config.get("mode", "build") != "build":
                        raise ValueError("switch to build mode before exporting a model")
                    current_input = self._active_input_path()
                    if current_input is None or (feature.key != "wyckoff" and not current_input.is_file()):
                        raise ValueError("choose an input structure before exporting")
                    if not self.active_feature.available:
                        raise ValueError("the selected feature is not available")
                    if feature.key == "wyckoff":
                        self.wyckoff_panel.backend.parameters(config)
                        if not (
                            cached is not None
                            and self.rendered_feature_key == "wyckoff"
                            and self.rendered_config == config
                            and self.rendered_input_path == current_input
                        ):
                            raise ValueError("click Generate crystal for the current Wyckoff parameters before exporting")
                    elif feature.key == "interface":
                        self.interface_panel.backend.parameters(config)
                    else:
                        _omnigb_parameters(config)
                    input_path = current_input
                    if not (
                        self.rendered_config == config
                        and self.rendered_input_path == input_path
                        and self.rendered_feature_key == self.active_feature.key
                    ):
                        cached = None
                model_gb = _load_model_module(feature.script_path or MODEL_SCRIPT)
                frame_count = (
                    len(cached.frames) if cached is not None
                    else math.prod(_load_model_module().cnid_grid(config))
                )
            except (OSError, RuntimeError, json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"Could not export current settings: {exc}")
                return
            dialog = StructureExportDialog(self, frame_index, frame_count)
            if dialog.exec() != QDialog.DialogCode.Accepted:
                return
            scope = str(dialog.scope_combo.currentData())
            fmt = str(dialog.format_combo.currentData())
            base_name = (
                input_path.stem if source is not None
                else f"wyckoff_SG{config['space_group']}" if feature.key == "wyckoff"
                else model_gb._default_output_name("build", config)
            )
            if scope == "all":
                default_name = f"{base_name}_all_frames.xyz"
            else:
                extension = {"poscar": ".vasp", "cif": ".cif", "extxyz": ".xyz"}[fmt]
                default_name = f"{base_name}_frame_{frame_index}{extension}"
            filters = {
                "poscar": "POSCAR (*.vasp POSCAR*)",
                "cif": "CIF (*.cif)",
                "extxyz": "Extended XYZ (*.xyz *.extxyz)",
            }
            selected, _ = QFileDialog.getSaveFileName(
                self, "Export structure / trajectory",
                str((self._dialog_directory(input_path.parent) if feature.key == "wyckoff"
                     else input_path.parent) / default_name), filters[fmt],
            )
            if not selected:
                return
            target = Path(selected).expanduser().resolve()
            if fmt == "extxyz" and target.suffix.lower() not in {".xyz", ".extxyz"}:
                target = target.with_name(target.name + ".xyz")
            if fmt == "cif" and target.suffix.lower() != ".cif":
                target = target.with_name(target.name + ".cif")
            self._remember_file_directory(target)
            request = StructureExportRequest(
                target, fmt, input_path, config, feature, scope, frame_index
            )
            if cached is not None:
                self._write_structure_export(cached, request)
                return
            self._export_request = request
            self._set_export_enabled(False)
            self.error_label.clear()
            self.render_status.setText("Building current CNID grid for export…")
            self._build_timer.stop()
            self._start_export_build()

        def _set_export_enabled(self, enabled: bool) -> None:
            """Keep toolbar and menu export actions in the same state."""

            self.export_action.setEnabled(enabled)
            self.export_menu_action.setEnabled(enabled)

        def _start_export_build(self) -> None:
            """Build the saved export request after other backend work finishes."""

            request = self._export_request
            if request is None or self._export_build_active or self._export_write_active:
                return
            if (
                self.worker is not None and self.worker.isRunning()
            ) or (
                self.scan_worker is not None and self.scan_worker.isRunning()
            ):
                self._export_timer.start()
                return
            self._export_timer.stop()
            self._build_timer.stop()
            self.build_serial += 1
            self._export_build_active = True
            self.render_status.setText("Building current parameters for export…")
            self.worker = BuildWorker(
                request.input_path, request.config, request.feature
            )
            self.worker.completed.connect(
                lambda structure, snapshot=request: self._export_ready(
                    structure, snapshot
                )
            )
            self.worker.failed.connect(
                lambda error, snapshot=request: self._export_failed(error, snapshot)
            )
            serial = self.build_serial
            self.worker.progress.connect(
                lambda completed, total, token=serial:
                self.build_progress(token, completed, total)
            )
            self.worker.start()

        def _export_ready(
            self, structure: Any, request: StructureExportRequest
        ) -> None:
            """Save the requested model, regardless of subsequent control edits."""

            if self._export_request is not request:
                return
            self._export_request = None
            self._export_build_active = False
            self._export_timer.stop()
            if self._request_matches_editor(
                request.config, request.input_path, request.feature.key
            ):
                # The exact export model also becomes the current preview.
                # A duplicate refresh may have changed the serial without
                # changing any settings, so accept this matching snapshot.
                self.preview_ready(
                    self.build_serial,
                    structure,
                    request.config,
                    request.input_path,
                    request.feature.key,
                )
            else:
                self.render_status.setText("Export model built from the saved parameters")
            if self.trajectory is structure and request.frame_index < len(structure.frames):
                self.frame_changed(request.frame_index)
            self._set_export_enabled(self.last_structure is not None)
            self._write_structure_export(structure, request)

        def _export_failed(
            self, error: str, request: StructureExportRequest
        ) -> None:
            """Report a failed current build without exporting an older model."""

            if self._export_request is not request:
                return
            self._export_request = None
            self._export_build_active = False
            self._export_timer.stop()
            self._set_export_enabled(self.last_structure is not None)
            self.render_status.setText("Export build failed")
            self.error_label.setText(
                f"Could not build current parameters; no structure was exported: {error}"
            )
            self._resume_preview_after_export(request)

        def _resume_preview_after_export(self, request: StructureExportRequest) -> None:
            """Resume live updates only if the user edited the export snapshot."""

            if self._closing:
                self._build_timer.stop()
                return
            if not self._request_matches_editor(
                request.config, request.input_path, request.feature.key
            ) and self._active_auto_refresh():
                self.queue_build()
            else:
                self._build_timer.stop()

        def _write_structure_export(
            self, structure: StructureTrajectory, request: StructureExportRequest
        ) -> None:
            """Write the captured frame or all frames without losing cell/PBC."""

            self._export_request = request
            self._export_write_active = True
            self._set_export_enabled(False)
            self.export_worker = StructureExportWorker(structure, request, self)
            worker = self.export_worker
            worker.finished.connect(lambda active=worker: self._finish_export_worker(active))
            self.export_worker.completed.connect(
                lambda snapshot=request, count=len(structure.frames):
                self._export_written(snapshot, count)
            )
            self.export_worker.failed.connect(
                lambda error, snapshot=request: self._export_write_failed(error, snapshot)
            )
            self.render_status.setText("Writing structure / trajectory…")
            self.export_worker.start()

        def _finish_export_worker(self, worker: StructureExportWorker) -> None:
            """Retain the writer until its native thread has finished."""

            if self.export_worker is worker:
                self.export_worker = None
            worker.deleteLater()

        def _export_written(self, request: StructureExportRequest, count: int) -> None:
            """Finish the captured export and resume updates after file writing."""

            if self._export_request is not request:
                return
            self._export_request = None
            self._export_write_active = False
            self._set_export_enabled(self.trajectory is not None)
            selection = (
                f"all {count} frames"
                if request.scope == "all" else f"frame {request.frame_index}"
            )
            self.error_label.setText(f"Exported {selection}: {request.target}")
            if self.trajectory is not None and not self._closing:
                self.frame_changed(self.current_frame_index)
            if self.trajectory_source is None:
                self._resume_preview_after_export(request)

        def _export_write_failed(self, error: str, request: StructureExportRequest) -> None:
            """Report an export failure; an existing destination stays intact."""

            if self._export_request is not request:
                return
            self._export_request = None
            self._export_write_active = False
            self._set_export_enabled(self.trajectory is not None)
            self.error_label.setText(f"Could not export structure / trajectory: {error}")
            if self.trajectory is not None and not self._closing:
                self.frame_changed(self.current_frame_index)
            if self.trajectory_source is None:
                self._resume_preview_after_export(request)

        def _read_editor_config(self) -> dict[str, Any]:
            """Return the current JSON editor contents as a configuration."""

            if self.active_feature.key == "wyckoff":
                return self.wyckoff_panel.config()
            if self.active_feature.key == "interface":
                return _viewer_config(self.interface_panel.config())
            self.rotation_axis_editor.commit()
            self.gb_plane_editor.commit()
            # Toolbar/menu actions can run before a focused spinbox commits
            # its text. Finish those edits before reading the live JSON.
            for spin_type in (QDoubleSpinBox, QSpinBox):
                for spin in self.omnigb_panel.findChildren(spin_type):
                    spin.interpretText()
            config = json.loads(self.json_editor.toPlainText())
            if not isinstance(config, dict):
                raise ValueError("gbin.json must contain a JSON object")
            return _viewer_config(config)

        def save_config(self) -> None:
            if self.active_feature.key == "wyckoff":
                self.wyckoff_panel.save_config()
                return
            if self.active_feature.key == "interface":
                self.interface_panel.save_config()
                return
            if self.config_path is None:
                directory = self._dialog_directory(
                    self.input_path.parent if self.input_path else REPOSITORY_ROOT
                )
                selected, _ = QFileDialog.getSaveFileName(
                    self,
                    "Save gbin.json",
                    str(directory / "gbin.json"),
                    "JSON files (*.json)",
                )
                if not selected:
                    return
                self.config_path = Path(selected).resolve()
                self._remember_file_directory(self.config_path)
                self.config_label.setText(str(self.config_path))
            try:
                config = self._read_editor_config()
                _write_config(self.config_path, config)
                self._remember_file_directory(self.config_path)
                self.error_label.setText(f"Saved: {self.config_path}")
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"Could not save gbin.json: {exc}")

        def save_config_as(self) -> None:
            """Save the live JSON editor state to a separate snapshot file."""

            if self.active_feature.key == "wyckoff":
                self.wyckoff_panel.save_config()
                return
            if self.active_feature.key == "interface":
                self.interface_panel.save_config()
                return
            try:
                config = self._read_editor_config()
            except (json.JSONDecodeError, TypeError, ValueError) as exc:
                self.error_label.setText(f"Could not save current JSON: {exc}")
                return
            if self.config_path is not None:
                fallback_directory = self.config_path.parent
            elif self.input_path is not None:
                fallback_directory = self.input_path.parent
            else:
                fallback_directory = REPOSITORY_ROOT
            default_path = self._dialog_directory(fallback_directory) / "gbin.json"
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Save current JSON as",
                str(default_path),
                "JSON files (*.json);;All files (*)",
            )
            if not selected:
                return
            target = Path(selected).expanduser().resolve()
            self._remember_file_directory(target)
            try:
                _write_config(target, config)
            except (OSError, TypeError, ValueError) as exc:
                self.error_label.setText(
                    f"Could not save current JSON: {exc}"
                )
                return
            self.error_label.setText(f"Saved current JSON: {target}")

        def closeEvent(self, event: Any) -> None:
            running_worker = next(
                (
                    worker for worker in (self.worker, self.scan_worker, self.export_worker)
                    if worker is not None and worker.isRunning()
                ),
                None,
            )
            self._closing = True
            if hasattr(self, "interface_panel"):
                self.interface_panel._timer.stop()
            self._export_timer.stop()
            if not self._export_build_active and not self._export_write_active:
                self._export_request = None
            if running_worker is not None:
                # InterfaceMaster cannot be interrupted halfway through a build.
                # Defer destruction so Qt never deletes a running QThread.
                self.auto_refresh = False
                self._build_timer.stop()
                self._json_timer.stop()
                self.render_status.setText("Finishing the build before closing…")
                if isinstance(running_worker, BuildWorker):
                    running_worker.requestInterruption()
                running_worker.finished.connect(self.close)
                event.ignore()
                return
            self.viewport.close()
            event.accept()

    return QApplication, GrainBoundaryViewer


def _default_paths() -> tuple[Path | None, Path | None]:
    """Find the example files without requiring a working directory."""

    if getattr(sys, "frozen", False):
        # Double-click opens the file picker UI instead of starting a bundled
        # example calculation before the user has chosen their own files.
        return None, None
    directory = REPOSITORY_ROOT / "example"
    input_path = directory / "ZrO2.cif"
    config_path = directory / "gbin.json"
    return (
        input_path if input_path.is_file() else None,
        config_path if config_path.is_file() else None,
    )


def _verify_window(
    application: Any, window: Any, directory: Path, expect_preview: bool
) -> None:
    """Verify the actual native window and save its rendered contents."""

    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QToolBar

    directory.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + 120.0
    timer = QTimer(window)
    timer.setInterval(500)

    def check() -> None:
        handle = window.windowHandle()
        exposed = handle is not None and handle.isExposed()
        worker_running = window.worker is not None and window.worker.isRunning()
        ready = not expect_preview or window.last_structure is not None
        failed = bool(window.error_label.text()) and not worker_running
        if time.monotonic() < deadline and not failed:
            if not exposed or not ready or worker_running:
                return
        timer.stop()
        actions = [
            action.text()
            for toolbar in window.findChildren(QToolBar)
            for action in toolbar.actions()
            if not action.isSeparator()
        ]
        result = {
            "success": bool(exposed and ready and not failed),
            "window_visible": window.isVisible(),
            "window_exposed": bool(exposed),
            "toolbar_actions": actions,
            "atoms": len(window.last_structure) if ready and expect_preview else 0,
            "status": window.render_status.text(),
            "error": window.error_label.text(),
        }
        if expect_preview and ready:
            try:
                camera_before = window.viewport.camera
                camera_state_before = (
                    tuple(camera_before.position),
                    tuple(camera_before.focal_point),
                    tuple(camera_before.up),
                    (float(camera_before.parallel_scale),),
                    (float(camera_before.view_angle),),
                )
                window._render_structure(window.last_structure)
                camera_after = window.viewport.camera
                camera_state_after = (
                    tuple(camera_after.position),
                    tuple(camera_after.focal_point),
                    tuple(camera_after.up),
                    (float(camera_after.parallel_scale),),
                    (float(camera_after.view_angle),),
                )
                result["camera_preserved"] = all(
                    abs(left - right) < 1e-6
                    for before, after in zip(camera_state_before, camera_state_after)
                    for left, right in zip(before, after)
                )
                export_path = directory / "feature_check.POSCAR"
                window.last_structure.to(fmt="poscar", filename=str(export_path))
                result["export_roundtrip"] = export_path.is_file()
                import ase
                import numpy as np

                trajectory = window.trajectory
                if trajectory is None:
                    raise RuntimeError("the preview did not return a trajectory")
                expected_count = math.prod(
                    _load_model_module().cnid_grid(window.rendered_config or {})
                )
                result["grid_frame_count"] = len(trajectory.frames) == expected_count
                result["frames"] = len(trajectory.frames)
                selected_index = len(trajectory.frames) - 1
                window.frame_slider.setValue(selected_index)
                result["frame_slider_sync"] = (
                    window.frame_spin.value() == selected_index
                    and window.current_frame_index == selected_index
                )
                window.frame_spin.setValue(0)
                result["frame_editor_sync"] = (
                    window.frame_slider.value() == 0
                    and window.current_frame_index == 0
                )
                all_request = StructureExportRequest(
                    directory / "all_frames.xyz", "extxyz",
                    window.rendered_input_path, window.rendered_config or {},
                    window.active_feature, "all", selected_index,
                )
                _export_trajectory(trajectory, all_request)
                restored = _read_trajectory(all_request.target)
                result["trajectory_roundtrip"] = (
                    len(restored.frames) == len(trajectory.frames)
                    and all(
                        [str(site.specie) for site in source]
                        == [str(site.specie) for site in reread]
                        and np.allclose(source.cart_coords, reread.cart_coords, atol=1e-7)
                        and np.allclose(source.lattice.matrix, reread.lattice.matrix, atol=1e-7)
                        and source.lattice.pbc == reread.lattice.pbc
                        for source, reread in zip(trajectory.frames, restored.frames)
                    )
                )
                current_request = StructureExportRequest(
                    directory / "current_frame.xyz", "extxyz",
                    window.rendered_input_path, window.rendered_config or {},
                    window.active_feature, "current", selected_index,
                )
                _export_trajectory(trajectory, current_request)
                current = _read_trajectory(current_request.target)
                result["current_frame_roundtrip"] = (
                    len(current.frames) == 1
                    and np.allclose(
                        current.frames[0].cart_coords,
                        trajectory.frames[selected_index].cart_coords, atol=1e-7,
                    )
                )
                result["ase_version"] = ase.__version__
                result["ase_module_path"] = str(ase.__file__)
                result["frozen"] = bool(getattr(sys, "frozen", False))
                if result["frozen"]:
                    result["ase_bundled"] = Path(ase.__file__).resolve().is_relative_to(
                        REPOSITORY_ROOT.resolve()
                    )
                species = [str(site.specie) for site in window.last_structure]
                centers, neighbours, *_ = window.last_structure.get_neighbor_list(3.0)
                pair_names = sorted(
                    {
                        "-".join(sorted((species[int(first)], species[int(second)])))
                        for first, second in zip(centers, neighbours)
                    }
                )
                result["pair_filter"] = False
                if pair_names:
                    import numpy as np
                    import pyvista as pv

                    baseline = _bond_lines(
                        np, pv, window.last_structure, 0.0, 3.0
                    )
                    filtered = _bond_lines(
                        np,
                        pv,
                        window.last_structure,
                        0.0,
                        3.0,
                        {pair_names[0]: {"enabled": False}},
                    )
                    result["pair_filter"] = filtered.n_cells < baseline.n_cells
            except Exception as feature_error:
                result["feature_error"] = str(feature_error)
        if expect_preview:
            checks = (
                "camera_preserved", "export_roundtrip", "pair_filter",
                "grid_frame_count", "frame_slider_sync", "frame_editor_sync",
                "trajectory_roundtrip", "current_frame_roundtrip",
            )
            result["success"] = result["success"] and all(
                result.get(name, False) for name in checks
            ) and "feature_error" not in result
            if getattr(sys, "frozen", False):
                result["success"] = result["success"] and result.get("ase_bundled", False)
        result["screenshot_saved"] = window.grab().save(
            str(directory / "window.png")
        )
        if expect_preview and ready:
            window.viewport.screenshot(str(directory / "crystal.png"))
        (directory / "verification.json").write_text(
            json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        application.setProperty("verification_failed", not result["success"])
        window.close()

    timer.timeout.connect(check)
    timer.start()


def main() -> int:
    """Launch the Windows desktop viewer."""

    parser = argparse.ArgumentParser(
        description="Interactive CCkit grain-boundary viewer for model-gb.py."
    )
    parser.add_argument(
        "structure",
        nargs="?",
        help="crystal structure file (CIF, POSCAR, or another pymatgen format)",
    )
    parser.add_argument(
        "config",
        nargs="?",
        help="gbin.json configuration file",
    )
    parser.add_argument("--verify", type=Path, help=argparse.SUPPRESS)
    arguments = parser.parse_args()

    # Keep the stream alive until Qt exits, including during VTK native calls.
    diagnostic_stream = _configure_diagnostics()
    diagnostic_stream.write("Loading Qt and renderer dependencies…\n")
    try:
        QApplication, viewer_type = create_application_classes()
    except RuntimeError as exc:
        _report_startup_error(exc)
        return 2

    application = QApplication(sys.argv)
    application.setApplicationName("CCkit Grain Boundary Viewer")
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QPixmap
    from PySide6.QtWidgets import QSplashScreen

    splash_image = QPixmap(460, 160)
    splash_image.fill(Qt.GlobalColor.darkBlue)
    splash = QSplashScreen(splash_image)
    splash.showMessage(
        "CCkit Grain Boundary Viewer\nLoading the renderer…",
        Qt.AlignmentFlag.AlignCenter,
        Qt.GlobalColor.white,
    )
    splash.show()
    application.processEvents()

    input_path, config_path = _default_paths()
    if arguments.structure:
        input_path = Path(arguments.structure).expanduser().resolve()
    if arguments.config:
        config_path = Path(arguments.config).expanduser().resolve()

    window = viewer_type(input_path, config_path)
    window.show()
    splash.finish(window)
    diagnostic_stream.write("Main window shown.\n")
    if arguments.verify is not None:
        _verify_window(
            application, window, arguments.verify, input_path is not None
        )
    exit_code = application.exec()
    diagnostic_stream.write(f"Application exited: {exit_code}\n")
    if application.property("verification_failed"):
        return 1
    return exit_code


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        _report_startup_error(error)
        raise SystemExit(2)
