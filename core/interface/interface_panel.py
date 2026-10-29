"""Independent engine selection, ZSL matching and atomic interface controls."""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Callable

from PySide6.QtCore import QPointF, QRectF, QThread, QTimer, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (
    QCheckBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLineEdit,
    QPlainTextEdit, QPushButton, QSizePolicy, QVBoxLayout, QWidget,
)

from viewer_widgets import OneDecimalDisplaySpinBox, WrappingLabel as QLabel


class MatchPreview(QWidget):
    """Compare reduced surface cells on one scale, with first vectors aligned."""

    def __init__(self) -> None:
        super().__init__()
        self.record: dict[str, Any] | None = None
        self.setMinimumHeight(170)
        self.setMaximumHeight(210)
        self.setMinimumWidth(120)
        self.setToolTip("Reduced surface cells before strain, with their first vectors aligned on a common length scale.")

    def set_record(self, record: dict[str, Any] | None) -> None:
        self.record = record
        self.update()

    def paintEvent(self, event: Any) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#f8fafc"))
        painter.setPen(QColor("#475569"))
        if self.record is None:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Search and select a matching cell")
            return
        polygons = []
        for side in ("substrate", "film"):
            a, b = self.record[f"{side}_vectors_A"]
            length_a = math.sqrt(sum(x * x for x in a))
            dot = sum(x * y for x, y in zip(a, b)) / length_a
            height = math.sqrt(max(0, sum(x * x for x in b) - dot * dot))
            polygons.append([(0.0, 0.0), (length_a, 0.0), (length_a + dot, height), (dot, height)])
        max_x = max(x for points in polygons for x, _ in points)
        min_x = min(x for points in polygons for x, _ in points)
        max_y = max(y for points in polygons for _, y in points)
        scale = min((self.width() - 48) / max(max_x - min_x, 1e-9), (self.height() - 65) / max(max_y, 1e-9))
        for points, color, style in zip(
            polygons, ("#d97706", "#2563eb"), (Qt.PenStyle.SolidLine, Qt.PenStyle.DashLine),
        ):
            painter.setPen(QPen(QColor(color), 2, style))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPolygon(QPolygonF([
                QPointF(24 + (x - min_x) * scale, self.height() - 42 - y * scale)
                for x, y in points
            ]))
        painter.setPen(QColor("#d97706"))
        painter.drawText(12, 18, "Substrate")
        painter.setPen(QColor("#2563eb"))
        painter.drawText(110, 18, "Film")
        painter.setPen(QColor("#475569"))
        painter.drawText(12, self.height() - 12, "Before strain · a aligned · equal scale")


class DichromaticPreview(MatchPreview):
    """Show lattice translations in the strained CSL and its DSCL cell."""

    def __init__(self) -> None:
        super().__init__()
        self.points: dict[str, list[list[float]]] = {}
        self.setMinimumHeight(240)
        self.setMaximumHeight(300)
        self.setToolTip("Ideal commensurate Bravais sites after film deformation. Chemical atomic positions are shown in the main 3D view after building. At most 250 sites per lattice are drawn.")

    def set_record(self, record: dict[str, Any] | None) -> None:
        from core.interface.commensurability import dichromatic_points

        self.points = dichromatic_points(record) if record is not None else {}
        super().set_record(record)

    def paintEvent(self, event: Any) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#f8fafc"))
        painter.setPen(QColor("#475569"))
        if self.record is None:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "CSL / DSCL dichromatic pattern")
            return
        a, b = (vector[:2] for vector in self.record["substrate_vectors_A"])
        da, db = self.record["dscl_vectors_A"]
        cell = [(0, 0), a, [a[0] + b[0], a[1] + b[1]], b]
        dscl = [(0, 0), da, [da[0] + db[0], da[1] + db[1]], db]
        if self.width() >= 420:
            main_rect = QRectF(12, 30, self.width() - 204, self.height() - 62)
            inset_rect = QRectF(self.width() - 180, 30, 168, self.height() - 62)
        else:
            main_rect = QRectF(12, 28, self.width() - 24, (self.height() - 60) * 0.5)
            inset_rect = QRectF(12, main_rect.bottom() + 8, self.width() - 24, (self.height() - 60) * 0.5)

        def mapper(polygon: Any, rectangle: QRectF) -> Callable[[Any], QPointF]:
            min_x, max_x = min(p[0] for p in polygon), max(p[0] for p in polygon)
            min_y, max_y = min(p[1] for p in polygon), max(p[1] for p in polygon)
            scale = min(rectangle.width() / max(max_x - min_x, 1e-9),
                        rectangle.height() / max(max_y - min_y, 1e-9))
            center = rectangle.center()
            return lambda xy: QPointF(center.x() + (xy[0] - (min_x + max_x) / 2) * scale,
                                     center.y() - (xy[1] - (min_y + max_y) / 2) * scale)

        point = mapper(cell + dscl, main_rect.adjusted(6, 6, -6, -6))

        for polygon, color in ((cell, "#64748b"), (dscl, "#9333ea")):
            painter.setPen(QPen(QColor(color), 1.5))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPolygon(QPolygonF([point(p) for p in polygon]))
        for side, color in (("substrate", "#d97706"), ("film", "#2563eb")):
            painter.setPen(QPen(QColor(color), 1.4))
            for xy in self.points[side]:
                p = point(xy)
                if side == "substrate":
                    painter.drawEllipse(p, 3.3, 3.3)
                else:
                    painter.drawLine(p + QPointF(-2, -2), p + QPointF(2, 2))
                    painter.drawLine(p + QPointF(-2, 2), p + QPointF(2, -2))
        painter.setPen(QPen(QColor("#cbd5e1"), 1))
        painter.drawRect(inset_rect)
        painter.setPen(QColor("#9333ea"))
        painter.drawText(inset_rect.adjusted(4, 2, -4, -2), Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignHCenter,
                         "DSCL · separate scale")
        inset_point = mapper(dscl, inset_rect.adjusted(14, 27, -14, -37))
        painter.setPen(QPen(QColor("#9333ea"), 2))
        painter.drawPolygon(QPolygonF([inset_point(p) for p in dscl]))
        painter.setBrush(QColor("#9333ea"))
        painter.drawEllipse(inset_point((0, 0)), 2.5, 2.5)
        lengths = [math.hypot(*vector) for vector in (da, db)]
        painter.drawText(inset_rect.adjusted(4, 2, -4, -3), Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
                         f"|d₁| {lengths[0]:.3f} Å\n|d₂| {lengths[1]:.3f} Å")
        painter.setPen(QColor("#475569"))
        painter.drawText(10, 18, "CSL · substrate ○ / film ×")
        painter.drawText(10, self.height() - 8, "Ideal Bravais sites · up to 250 per lattice")


class MatchSearchWorker(QThread):
    """Keep native lattice searches off the Qt event loop."""

    completed = Signal(object)
    failed = Signal(str)

    def __init__(self, backend: Any, path: Path, config: dict[str, Any], parent: QWidget) -> None:
        super().__init__(parent)
        self.backend, self.path, self.config = backend, path, copy.deepcopy(config)

    def run(self) -> None:
        try:
            result = self.backend.search_matches(
                self.path, self.config, cancelled=self.isInterruptionRequested,
            )
        except Exception as exc:  # Report native-library errors at the worker boundary.
            self.failed.emit(str(exc))
        else:
            self.completed.emit(result)


class InterfacePanel(QWidget):
    """Own engine settings, candidates and JSON independently of other features."""

    changed = Signal()
    build_requested = Signal()
    failed = Signal(str)

    def __init__(
        self, backend: Any, index_editor: Any, integer_spin: Any,
        double_spin: Any, combo_box: Any, file_filter: str,
        dialog_directory: Callable[[Path], Path], remember_directory: Callable[[Path], None],
    ) -> None:
        super().__init__()
        self.backend = backend
        self.input_path: Path | None = None
        self.config_path: Path | None = None
        self._directory, self._remember, self._file_filter = dialog_directory, remember_directory, file_filter
        self._updating = True
        self._state = backend.default_config()
        self._result: dict[str, Any] | None = None
        self._last_report: dict[str, Any] | None = None
        self._result_stamp: str | None = None
        self._search_serial = 0
        self.search_worker: MatchSearchWorker | None = None
        self.numbers: dict[str, Any] = {}
        self.checks: dict[str, QCheckBox] = {}
        self.indices: dict[str, Any] = {}
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        engine_group = QGroupBox("Modeling engine")
        engine_form = QFormLayout(engine_group)
        engine_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.engine_combo = combo_box()
        self.engine_combo.setMinimumContentsLength(12)
        self.engine_combo.setSizeAdjustPolicy(self.engine_combo.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.engine_combo.addItem("pymatgen · Zur–McGill", "pymatgen")
        self.engine_combo.addItem("interface-master · CSL / DSCL", "interface-master")
        self.engine_combo.currentIndexChanged.connect(self.controls_changed)
        engine_form.addRow("Engine", self.engine_combo)
        layout.addWidget(engine_group)
        self.engine_note = QLabel()
        self.engine_note.setWordWrap(True)
        layout.addWidget(self.engine_note)
        self.content = QWidget()
        content = QVBoxLayout(self.content)
        content.setContentsMargins(0, 0, 0, 0)
        content.setSpacing(8)
        layout.addWidget(self.content)

        def group(title: str) -> QFormLayout:
            box = QGroupBox(title)
            form = QFormLayout(box)
            form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
            form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
            content.addWidget(box)
            return form

        def number(form: QFormLayout, key: str, label: str, lower: float, upper: float,
                   integer: bool = False, tooltip: str = "") -> None:
            display_only_precision = key in {"gap", "offset_a", "offset_b"}
            if display_only_precision:
                spin = OneDecimalDisplaySpinBox()
                tooltip = (
                    f"{tooltip}\nDisplayed to one decimal; calculations and JSON retain "
                    "the full entered value."
                ).strip()
            else:
                spin = integer_spin() if integer else double_spin()
            spin.setRange(int(lower), int(upper)) if integer else spin.setRange(lower, upper)
            if not integer:
                if not display_only_precision:
                    spin.setDecimals(6)
                spin.setSingleStep(0.01 if "tol" in key or "offset" in key
                                   or key in {"im_max_strain", "im_max_area_change"} else 0.5)
            spin.setKeyboardTracking(False)
            spin.setToolTip(tooltip)
            spin.valueChanged.connect(self.controls_changed)
            self.numbers[key] = spin
            form.addRow(label, spin)

        def check(form: QFormLayout, key: str, label: str, tooltip: str = "") -> None:
            widget = QCheckBox(label)
            widget.setToolTip(tooltip)
            widget.toggled.connect(self.controls_changed)
            self.checks[key] = widget
            form.addRow(widget)

        files = group("1 · Input structures and contact planes")
        self.substrate_label, self.film_label = QLabel("not selected"), QLabel("not selected")
        for side, label in (("substrate", self.substrate_label), ("film", self.film_label)):
            label.setWordWrap(True)
            label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            button = QPushButton(f"Choose {side}…")
            button.clicked.connect(lambda _checked=False, side=side: self.choose_structure(side))
            files.addRow(button, label)
            editor = index_editor("plane")
            editor.changed.connect(self.controls_changed)
            editor.invalid.connect(self.failed.emit)
            self.indices[f"{side}_plane"] = editor
            files.addRow(f"{side.title()} plane", editor)
        hint = QLabel(
            "Indices use each input cell's basis. "
            "Miller–Bravais indices require a standard hexagonal cell."
        )
        hint.setWordWrap(True)
        files.addRow(hint)

        search = group("2 · Surface lattice matching")
        number(search, "max_area", "Maximum area (Å²)", 1, 100000,
               tooltip="Both unstrained surface supercell areas must stay within this bound.")
        self.zsl_controls = QWidget()
        zsl = QFormLayout(self.zsl_controls)
        zsl.setContentsMargins(0, 0, 0, 0)
        zsl.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        search.addRow(self.zsl_controls)
        for key, label, tip in (
            ("max_area_ratio_tol", "Area-ratio filter", "Native max_area_ratio_tol: filters candidate area multiples; 0.09 is not a final 9% strain limit."),
            ("max_length_tol", "Relative length tolerance", "Native max_length_tol: |substrate length / film length − 1|; 0.03 = 3%."),
            ("max_angle_tol", "Relative angle tolerance", "Native max_angle_tol: |substrate angle / film angle − 1|; 0.01 = 1%, not 0.01 degrees."),
        ):
            number(zsl, key, label, 0, 1, tooltip=tip)
        check(zsl, "bidirectional", "Accept matching in either comparison direction",
              "Uses pymatgen's OR of film→substrate and substrate→film tolerance checks.")
        check(zsl, "lowest", "Retain only the first native match",
              "Native lowest=True stops at the first match in its search order; this is not an energy ranking.")
        self.im_controls = QWidget()
        im_form = QFormLayout(self.im_controls)
        im_form.setContentsMargins(0, 0, 0, 0)
        im_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        for key, label, lower, upper, integer, tip in (
            ("im_angle_start", "Rotation start (°)", -360, 360, False, "Normals are aligned; 0° aligns the primitive surface a directions."),
            ("im_angle_stop", "Rotation stop (°)", -360, 360, False, "Both endpoints are included. Choose the range for your material symmetry; there is no automatic symmetry reduction."),
            ("im_angle_step", "Rotation step (°)", 0.001, 360, False, "Explicit finite angular mesh; smaller steps explore orientations more finely."),
            ("im_rotation_tol", "Maximum rotation correction (°)", 0, 10, False, "Rational approximation can add a small polar rotation. Final angles must also lie in the requested range."),
            ("im_max_denominator", "Maximum rational denominator", 1, 100, True, "Search round(q U) / q for each q up to this bound; q is not the number of atoms or a CSL index."),
            ("im_max_strain", "Maximum principal strain", 0, 0.5, False, "Maximum absolute in-plane principal stretch minus one. 0.03 means 3%."),
            ("im_max_area_change", "Maximum film area change", 0, 1, False, "Require |det(D) − 1| below this threshold."),
            ("im_max_sigma", "Maximum commensurability σ", 1, 10000, False, "Inoue Eq. 27–28: sqrt(I_substrate × I_film) = sqrt(A_CSL / A_DSCL), for the strained approximant."),
            ("im_max_index", "Maximum index on either side", 1, 10000, True, "Limit the number of primitive surface cells in the minimal CSL on either side."),
        ):
            number(im_form, key, label, lower, upper, integer, tip)
        self.rank_combo = combo_box()
        self.rank_combo.addItem("Commensurability σ, then strain", "sigma")
        self.rank_combo.addItem("Principal strain, then σ", "strain")
        self.rank_combo.addItem("CSL area, then strain", "area")
        self.rank_combo.setMinimumContentsLength(12)
        self.rank_combo.setSizeAdjustPolicy(self.rank_combo.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.rank_combo.currentIndexChanged.connect(self.controls_changed)
        im_form.addRow("Rank candidates by", self.rank_combo)
        search.addRow(self.im_controls)
        number(search, "max_matches", "Retained candidate limit", 1, 10000, True,
               "InterfaceMaster scans the whole requested mesh and retains the best candidates. pymatgen stops retaining native matches at this limit.")
        for side in ("substrate", "film"):
            number(search, f"{side}_termination_ftol", f"{side.title()} plane tolerance (Å)", 0.000001, 5,
                   tooltip="Native termination_ftol: group nearby atoms into terminating planes.")
        check(search, "label_index", "Keep separately indexed termination labels")
        check(search, "filter_out_sym_slabs", "Filter symmetry-equivalent slabs")
        search_buttons = QHBoxLayout()
        self.search_button = QPushButton("Search matching cells")
        self.search_button.clicked.connect(self.search)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel_search)
        self.cancel_button.setEnabled(False)
        search_buttons.addWidget(self.search_button)
        search_buttons.addWidget(self.cancel_button)
        search.addRow(search_buttons)
        self.search_status = QLabel("Choose both structures, then search.")
        self.search_status.setWordWrap(True)
        search.addRow(self.search_status)

        result_form = group("3 · Compare candidates and choose terminations")
        self.match_combo = combo_box()
        self.match_combo.setMinimumContentsLength(12)
        self.match_combo.setSizeAdjustPolicy(self.match_combo.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.match_combo.currentIndexChanged.connect(self.controls_changed)
        result_form.addRow("Matching cell", self.match_combo)
        self.match_info = QLabel("No search results")
        self.match_info.setWordWrap(True)
        self.match_info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        result_form.addRow(self.match_info)
        self.preview = MatchPreview()
        result_form.addRow(self.preview)
        self.dichromatic = DichromaticPreview()
        result_form.addRow(self.dichromatic)
        self.period_info = QPlainTextEdit()
        self.period_info.setReadOnly(True)
        self.period_info.setMinimumHeight(130)
        self.period_info.setMaximumHeight(190)
        result_form.addRow(self.period_info)
        self.matrix_info = QPlainTextEdit()
        self.matrix_info.setReadOnly(True)
        self.matrix_info.setMinimumHeight(125)
        self.matrix_info.setMaximumHeight(175)
        result_form.addRow(self.matrix_info)
        self.termination_combo = combo_box()
        self.termination_combo.setMinimumContentsLength(10)
        self.termination_combo.setSizeAdjustPolicy(self.termination_combo.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.termination_combo.currentIndexChanged.connect(self.controls_changed)
        result_form.addRow("Film / substrate termination", self.termination_combo)
        self.report_button = QPushButton("Save matching / model report…")
        self.report_button.clicked.connect(self.save_report)
        result_form.addRow(self.report_button)

        geometry = group("4 · Build the atomic interface")
        self.system_edit = QLineEdit()
        self.system_edit.setMaxLength(200)
        self.system_edit.setPlaceholderText("e.g. hZrO-Zr · empty: substrate-film formulas")
        self.system_edit.setToolTip("System name used in each exported XYZ frame's config_type.")
        self.system_edit.textEdited.connect(self.controls_changed)
        geometry.addRow("System label", self.system_edit)
        self.thickness_units = combo_box()
        self.thickness_units.addItem("Oriented unit planes", True)
        self.thickness_units.addItem("Minimum thickness (Å)", False)
        self.thickness_units.currentIndexChanged.connect(self.controls_changed)
        geometry.addRow("Thickness units", self.thickness_units)
        for side in ("substrate", "film"):
            number(geometry, f"{side}_thickness", f"{side.title()} thickness", 0.000001, 10000)
        number(geometry, "gap", "Interface gap (Å)", 0, 10000)
        number(geometry, "vacuum_over_film", "External vacuum (Å)", 0, 10000)

        shift_group = QGroupBox("In-plane shift")
        shift_form = QFormLayout(shift_group)
        shift_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.shift_mode = combo_box()
        self.shift_mode.setMinimumContentsLength(12)
        self.shift_mode.setSizeAdjustPolicy(self.shift_mode.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.shift_mode.addItem("Custom shift · single structure", "custom")
        self.shift_mode.addItem("Grid sampling · multiple structures", "sampling")
        self.shift_mode.currentIndexChanged.connect(self.controls_changed)
        shift_form.addRow("Shift mode", self.shift_mode)
        self.custom_shift_panel = QWidget()
        custom_form = QFormLayout(self.custom_shift_panel)
        custom_form.setContentsMargins(0, 0, 0, 0)
        custom_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.custom_basis_panel = QWidget()
        custom_basis_form = QFormLayout(self.custom_basis_panel)
        custom_basis_form.setContentsMargins(0, 0, 0, 0)
        custom_basis_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.custom_basis_combo = combo_box()
        self.custom_basis_combo.addItem("DSCL cell", "dscl")
        self.custom_basis_combo.addItem("Matched a/b supercell", "supercell")
        self.custom_basis_combo.currentIndexChanged.connect(self.controls_changed)
        custom_basis_form.addRow("Shift basis", self.custom_basis_combo)
        custom_form.addRow(self.custom_basis_panel)
        for key, label in (("offset_a", "Coefficient d1"), ("offset_b", "Coefficient d2")):
            number(custom_form, key, label, -100, 100,
                   tooltip="Dimensionless coefficient of the selected first/second basis vector; negative values are allowed.")
        shift_form.addRow(self.custom_shift_panel)
        self.sample_basis_panel = QWidget()
        basis_form = QFormLayout(self.sample_basis_panel)
        basis_form.setContentsMargins(0, 0, 0, 0)
        basis_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.shift_basis_combo = combo_box()
        self.shift_basis_combo.addItem("DSCL cell", "dscl")
        self.shift_basis_combo.addItem("Matched a/b supercell", "supercell")
        self.shift_basis_combo.currentIndexChanged.connect(self.controls_changed)
        basis_form.addRow("Sampling cell", self.shift_basis_combo)
        shift_form.addRow(self.sample_basis_panel)
        self.sample_shift_panel = QWidget()
        sample_form = QFormLayout(self.sample_shift_panel)
        sample_form.setContentsMargins(0, 0, 0, 0)
        sample_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        for key, label in (("offset_grid_a", "Shift samples along a"), ("offset_grid_b", "Shift samples along b")):
            number(sample_form, key, label, 1, 100, True)
        shift_form.addRow(self.sample_shift_panel)
        self.shift_hint = QLabel()
        self.shift_hint.setWordWrap(True)
        shift_form.addRow(self.shift_hint)
        geometry.addRow(shift_group)
        self.geometry_note = QLabel()
        self.geometry_note.setWordWrap(True)
        geometry.addRow(self.geometry_note)
        self.registry_controls = QGroupBox("Atomic contact geometry")
        registry = QFormLayout(self.registry_controls)
        registry.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        number(registry, "im_contact_depth", "Contact-layer depth (Å)", 0, 100,
               tooltip="Measure projected registry using atoms within this distance of each facing surface.")
        number(registry, "im_registry_tolerance", "Projected coincidence tolerance (Å)", 0, 10,
               tooltip="Fraction of film contact-layer atoms whose nearest substrate projection lies within this radius; ignores species and is not an energy.")
        number(registry, "im_min_distance", "Close-contact distance (Å)", 0, 10,
               tooltip="Report film atoms closer than this distance to any substrate atom; no atoms are deleted.")
        number(registry, "im_max_atoms", "Maximum atoms per model", 1, 200000, True)
        geometry.addRow(self.registry_controls)
        self.auto_check = QCheckBox("Auto rebuild selected model after geometry changes")
        self.auto_check.toggled.connect(self.controls_changed)
        geometry.addRow(self.auto_check)
        self.build_button = QPushButton("Build selected interface")
        self.build_button.clicked.connect(self.request_build)
        geometry.addRow(self.build_button)
        self.summary = QLabel("No atomic interface built yet.")
        self.summary.setWordWrap(True)
        geometry.addRow(self.summary)
        self.frame_summary = QLabel()
        self.frame_summary.setWordWrap(True)
        self.frame_summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        geometry.addRow(self.frame_summary)

        config_group = QGroupBox("interface.json · advanced settings")
        config_group.setCheckable(True)
        config_group.setChecked(False)
        config_layout = QVBoxLayout(config_group)
        config_body = QWidget()
        body_layout = QVBoxLayout(config_body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        self.editor = QPlainTextEdit()
        self.editor.setMinimumHeight(220)
        body_layout.addWidget(self.editor)
        buttons = QHBoxLayout()
        for label, callback in (("Apply JSON", self.apply_json), ("Load JSON…", self.load_config), ("Save JSON…", self.save_config)):
            button = QPushButton(label)
            button.clicked.connect(callback)
            buttons.addWidget(button)
        body_layout.addLayout(buttons)
        config_layout.addWidget(config_body)
        config_group.toggled.connect(config_body.setVisible)
        config_body.setVisible(False)
        layout.addWidget(config_group)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(700)
        self._timer.timeout.connect(self.apply_json)
        self.editor.textChanged.connect(self.schedule_json)
        self._set_config(self._state)
        self.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)

    def _fingerprint(self, data: dict[str, Any]) -> str:
        return json.dumps([str(self.input_path), {key: data[key] for key in self.backend.SEARCH_KEYS}], sort_keys=True)

    def _file_stamp(self) -> str:
        """Invalidate candidates if either source file has changed externally."""

        if self.input_path is None:
            return ""
        film = Path(self._state["film_structure"]).expanduser()
        if not film.is_absolute():
            film = self.input_path.parent / film
        try:
            return repr([(str(p.resolve()), p.stat().st_mtime_ns, p.stat().st_size) for p in (self.input_path, film)])
        except OSError:
            return ""

    def ready_to_build(self) -> bool:
        return (
            self._result is not None
            and self._result_stamp == self._file_stamp()
            and 0 <= self._state["match_index"] < len(self._result["matches"])
            and self._state["termination"] in self._result["terminations"]
        )

    def _invalidate(self) -> None:
        self._search_serial += 1
        if self.search_worker is not None and self.search_worker.isRunning():
            self.search_worker.requestInterruption()
        self._result, self._result_stamp = None, None
        self._last_report = None
        previous = self._updating
        self._updating = True
        self.match_combo.clear()
        self.termination_combo.clear()
        self._updating = previous
        self.preview.set_record(None)
        self.dichromatic.set_record(None)
        self.period_info.clear()
        self.match_info.setText("Search again after changing structures, planes or matching parameters.")
        self.matrix_info.clear()
        self.search_status.setText("Search required for the current inputs.")
        self.summary.setText("The current settings have not been built.")
        self.frame_summary.clear()
        self._update_buttons()

    def _update_buttons(self) -> None:
        running = self.search_worker is not None and self.search_worker.isRunning()
        self.search_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)
        self.build_button.setEnabled(self.ready_to_build() and not running)
        self.report_button.setEnabled(self._result is not None and not running)
        im_engine = self._state["engine"] == "interface-master"
        self.zsl_controls.setVisible(not im_engine)
        self.im_controls.setVisible(im_engine)
        self.registry_controls.setVisible(im_engine)
        self.frame_summary.setVisible(True)
        self.dichromatic.setVisible(im_engine)
        self.period_info.setVisible(im_engine)
        self.engine_note.setText(
            "Approximate CSL / DSCL · rational matching periods inspired by Inoue et al. (2026). "
            "Compare finite rotation searches and build explicit atomic interfaces."
            if im_engine else "Zur–McGill (1984) surface superlattice matching with pymatgen."
        )
        self.geometry_note.setText(
            "The substrate is fixed; a bounded film deformation closes the periodic cell. "
            "Atomic registry and close contacts are geometric diagnostics. Relaxation and energy evaluation are separate steps. "
            "DSCL equivalence refers to ideal lattice translations; chemical terminations may have additional distinctions."
            if im_engine else "The substrate is the reference; the film is strained into coherency.\n"
            "Shift samples cover the matched a/b supercell, not a DSCL cell."
        )
        self._update_shift_controls()
        self.updateGeometry()

    def _update_shift_controls(self) -> None:
        custom = self._state["shift_mode"] == "custom"
        self.custom_shift_panel.setVisible(custom)
        self.sample_shift_panel.setVisible(not custom)
        im_engine = self._state["engine"] == "interface-master"
        self.sample_basis_panel.setVisible(im_engine and not custom)
        self.custom_basis_panel.setVisible(im_engine and custom)
        if custom:
            basis = self.backend.active_shift_basis(self._state)
            vectors = "DSCL vectors D1 / D2" if basis == "dscl" else "matched supercell vectors a / b"
            self.shift_hint.setText(
                f"Build one structure: translation = d1 × first vector + d2 × second vector, using {vectors}."
            )
        else:
            n1, n2 = self._state["offset_grid"]
            basis = ("DSCL cell" if im_engine and self._state["im_shift_basis"] == "dscl"
                     else "matched a/b supercell")
            self.shift_hint.setText(
                f"{n1} × {n2} = {n1 * n2} structures. Sampling starts at (0, 0); "
                f"positions are uniformly spaced over the {basis}."
            )

    def _store(self, data: dict[str, Any]) -> None:
        old_key = self._fingerprint(self._state)
        self._state = self.backend.parameters(data)
        if old_key != self._fingerprint(self._state):
            self._invalidate()
        previous = self._updating
        self._updating = True
        self.editor.setPlainText(json.dumps(self._state, indent=2, ensure_ascii=False))
        self._updating = previous
        self._update_buttons()
        self._show_match()
        self._timer.stop()

    def _thickness_options(self) -> None:
        layers = bool(self.thickness_units.currentData())
        for side in ("substrate", "film"):
            spin = self.numbers[f"{side}_thickness"]
            spin.setDecimals(0 if layers else 4)
            spin.setMinimum(1 if layers else 0.0001)
            spin.setSingleStep(1 if layers else 0.5)

    def controls_changed(self, *_args: Any) -> None:
        if self._updating or not hasattr(self, "editor"):
            return
        self._updating = True
        try:
            self._thickness_options()
            data = copy.deepcopy(self._state)
            data["engine"] = self.engine_combo.currentData()
            data["im_rank_by"] = self.rank_combo.currentData()
            data["im_shift_basis"] = self.shift_basis_combo.currentData()
            data["custom_shift_basis"] = self.custom_basis_combo.currentData()
            data["system"] = self.system_edit.text()
            for key, editor in self.indices.items():
                data[key], data[f"{key}_notation"] = editor.indices(), editor.notation()
            for key, spin in self.numbers.items():
                if not key.startswith("offset_grid_"):
                    data[key] = spin.value()
            data["offset_grid"] = [self.numbers["offset_grid_a"].value(), self.numbers["offset_grid_b"].value()]
            data.update({key: widget.isChecked() for key, widget in self.checks.items()})
            data["in_layers"] = bool(self.thickness_units.currentData())
            data["shift_mode"] = self.shift_mode.currentData()
            if self._result is not None:
                data["match_index"] = self.match_combo.currentData()
                data["termination"] = self.termination_combo.currentData()
            self._store(data)
        except (TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
        finally:
            self._updating = False
        self.changed.emit()

    def _set_config(self, config: dict[str, Any]) -> None:
        values = self.backend.parameters(config)
        if values.get("substrate_structure"):
            source = Path(values["substrate_structure"]).expanduser()
            if not source.is_absolute():
                directory = self.config_path.parent if self.config_path else (
                    self.input_path.parent if self.input_path else Path.cwd()
                )
                source = directory / source
            source = source.resolve()
            values["substrate_structure"] = str(source)
            if self.input_path != source:
                self.input_path = source
                self.substrate_label.setText(str(source))
                self._invalidate()
        self._updating = True
        try:
            self.engine_combo.setCurrentIndex(self.engine_combo.findData(values["engine"]))
            self.rank_combo.setCurrentIndex(self.rank_combo.findData(values["im_rank_by"]))
            self.shift_basis_combo.setCurrentIndex(self.shift_basis_combo.findData(values["im_shift_basis"]))
            self.custom_basis_combo.setCurrentIndex(self.custom_basis_combo.findData(values["custom_shift_basis"]))
            self.system_edit.setText(values["system"])
            self.shift_mode.setCurrentIndex(self.shift_mode.findData(values["shift_mode"]))
            for key, editor in self.indices.items():
                editor.set_indices(values[key], values[f"{key}_notation"])
            self.thickness_units.setCurrentIndex(self.thickness_units.findData(values["in_layers"]))
            self._thickness_options()
            for key, spin in self.numbers.items():
                value = values["offset_grid"][0 if key == "offset_grid_a" else 1] if key.startswith("offset_grid_") else values[key]
                spin.setValue(value)
            for key, check in self.checks.items():
                check.setChecked(values[key])
            self.film_label.setText(values["film_structure"] or "not selected")
            if self._result is not None:
                self.match_combo.setCurrentIndex(self.match_combo.findData(values["match_index"]))
                self.termination_combo.setCurrentIndex(self.termination_combo.findData(values["termination"]))
            self._store(values)
        finally:
            self._updating = False

    def config(self) -> dict[str, Any]:
        for spin in self.numbers.values():
            spin.interpretText()
        for editor in self.indices.values():
            editor.commit()
        values = self.backend.parameters(json.loads(self.editor.toPlainText()))
        if values != self._state:
            self._set_config(values)
        return copy.deepcopy(self._state)

    def schedule_json(self) -> None:
        if not self._updating:
            self._timer.start()

    def apply_json(self, *_args: Any) -> None:
        try:
            self._set_config(self.config())
        except (TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.changed.emit()

    def search(self) -> None:
        if self.search_worker is not None and self.search_worker.isRunning():
            return
        try:
            self._set_config(self.config())
            if self.input_path is None:
                raise ValueError("Choose a substrate structure first")
            if not self._state["film_structure"]:
                raise ValueError("Choose a film structure first")
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self._invalidate()
        serial = self._search_serial
        stamp = self._file_stamp()
        worker = MatchSearchWorker(self.backend, self.input_path, self._state, self)
        self.search_worker = worker
        worker.completed.connect(lambda result: self._search_ready(serial, stamp, result))
        worker.failed.connect(lambda error: self._search_failed(serial, error))
        worker.finished.connect(self._search_finished)
        self.search_status.setText("Searching surface superlattices and terminations…")
        worker.start()
        self._update_buttons()

    def cancel_search(self) -> None:
        self._invalidate()
        self.search_status.setText("Search cancelled; waiting for the current native step to finish…")
        self.changed.emit()

    def _search_ready(self, serial: int, stamp: str, result: dict[str, Any]) -> None:
        if serial != self._search_serial or stamp != self._file_stamp():
            return
        self._result, self._result_stamp = result, stamp
        self._updating = True
        for row in result["matches"]:
            label = (
                f"#{row['index'] + 1} · θ {row['twist_angle_deg']:.1f}° · σ {row['sigma']:.2f} · "
                f"strain {max(abs(x) for x in row['principal_strains']):.2%}"
                if "sigma" in row else
                f"#{row['index'] + 1} · {row['substrate_area_A2']:.2f} Å² · strain {row['von_mises_strain']:.4f}"
            )
            self.match_combo.addItem(
                label, row["index"],
            )
        index = self.match_combo.findData(self._state["match_index"])
        self.match_combo.setCurrentIndex(max(0, index))
        for pair in result["terminations"]:
            self.termination_combo.addItem(f"{pair[0]} / {pair[1]}", pair)
        index = self.termination_combo.findData(self._state["termination"])
        self.termination_combo.setCurrentIndex(max(0, index))
        self._updating = False
        count = len(result["matches"])
        suffix = " · candidate limit reached; increase limit for more results" if result["truncated"] else ""
        self.search_status.setText(f"{count} matching cells · {len(result['terminations'])} termination pairs{suffix}")
        if result.get("engine") == "interface-master":
            self.search_status.setText(
                f"{result['total_matches']} distinct rational candidates accepted over {result['trial_count']} trials; "
                f"{count} retained · {len(result['terminations'])} termination pairs. "
                "Search is bounded by the chosen angle mesh and rational denominators."
            )
            lines = ["1D period guides from primitive a/b lengths (not the rotated 2D solution)",
                     "M = film periods; N = substrate periods; CF = principal convergent"]
            for axis, rows in result["rational_periods"].items():
                lines.append(f"{axis} direction:")
                for period in rows:
                    parents = period["farey_parents"]
                    motif = (f" · {parents[0][0]}/{parents[0][1]} ⊕ {parents[1][0]}/{parents[1][1]}" if parents else "")
                    lines.append(f"  M:N {period['M']}:{period['N']} · Δ1 {period['delta_1']:+.3%} · Δ2 {period['delta_2']:.3%}"
                                 f"{' · CF' if period['convergent'] else ''}{motif}")
            self.period_info.setPlainText("\n".join(lines))
        self.controls_changed()

    def _search_failed(self, serial: int, error: str) -> None:
        if serial == self._search_serial:
            self.search_status.setText(error)
            self.failed.emit(error)

    def _search_finished(self) -> None:
        worker = self.search_worker
        self.search_worker = None
        if worker is not None:
            worker.deleteLater()
        self._update_buttons()

    def _show_match(self) -> None:
        if self._result is None or self.match_combo.currentIndex() < 0:
            return
        row = self._result["matches"][self.match_combo.currentData()]
        da, db = row["length_mismatch"]
        self.match_info.setText(
            f"Areas: substrate {row['substrate_area_A2']:.3f} / film {row['film_area_A2']:.3f} Å²\n"
            f"Length mismatch: a {da:+.3%} · b {db:+.3%}\n"
            f"Angle mismatch {row['angle_mismatch']:+.3%} ({row['angle_difference_deg']:+.4f}°)\n"
            f"Film von Mises strain: {row['von_mises_strain']:.6f}"
        )
        if "sigma" in row:
            principal = ", ".join(f"{value:+.4%}" for value in row["principal_strains"])
            self.match_info.setText(
                f"Relative rotation θ {row['twist_angle_deg']:.1f}° (mesh {row['sample_angle_deg']:.1f}°)\n"
                f"CSL indices: substrate {row['substrate_index']} · film {row['film_index']}\n"
                f"Commensurability σ {row['sigma']:.2f} · DSCL area {row['dscl_area_A2']:.6f} Å²\n"
                f"CSL area {row['substrate_area_A2']:.3f} Å² · principal film strains {principal}"
            )
            self.dichromatic.set_record(row)
        self.preview.set_record(row)
        lines = []
        for side in ("substrate", "film"):
            lines.append(f"{side.title()} supercell matrix (surface basis):")
            lines.extend("  " + str(vector) for vector in row[f"{side}_supercell_matrix"])
            lines.append("Supercell vectors (Å, common surface frame):" if "sigma" in row
                         else "Reduced supercell vectors (Å, bulk input frame):")
            lines.extend("  " + "  ".join(f"{x:.5f}" for x in vector) for vector in row[f"{side}_vectors_A"])
        if "sigma" in row:
            lines.append(f"Rational primitive-basis transform P/q; q = {row['rational_denominator']}:")
            lines.extend("  " + str(vector) for vector in row["rational_numerator"])
            lines.append("DSCL vectors (Å, common surface frame):")
            lines.extend("  " + str(vector) for vector in row["dscl_vectors_A"])
        self.matrix_info.setPlainText("\n".join(lines))

    def request_build(self) -> None:
        try:
            self._set_config(self.config())
            if not self.ready_to_build():
                raise ValueError("Search matching cells for the current inputs before building")
        except (TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.build_requested.emit()

    def show_report(self, details: dict[str, Any]) -> None:
        self._last_report = copy.deepcopy(details)
        shift_label = "Custom shift" if details["shift_mode"] == "custom" else "Grid sampling"
        self.summary.setText(
            f"{details['engine']} · match #{details['match_index'] + 1} · {details['atom_count']} atoms\n"
            f"Area {details['interface_area_A2']:.3f} Å² · film strain {details['von_mises_strain']:.6f}\n"
            f"{shift_label} · {details['frame_count']} structure(s)"
        )
        if details["engine"] == "interface-master":
            best = details["registry_best"]
            self.summary.setText(
                self.summary.text()
                + f"\nθ {details['twist_angle_deg']:.1f}° · σ {details['sigma']:.2f}"
                + f"\nHighest projected coincidence: frame {details['registry_best_frame']} · {best['registry_fraction']:.1%}"
                + f"\nMean lateral distance {best['mean_lateral_distance_A']:.4f} Å · minimum cross-interface distance {best['minimum_cross_distance_A']:.4f} Å"
                + f"\n{details['clashing_frames']} frame(s) have close contacts. Registry is species-independent geometry, not an energy ranking."
                + f"\n{details['dscl_method']}"
            )

    def save_report(self, *_args: Any) -> None:
        try:
            config = self.config()
            if self._result is None or self._result_stamp != self._file_stamp():
                raise ValueError("Search current inputs before saving a matching report")
            selected, _ = QFileDialog.getSaveFileName(
                self, "Save matching / model report",
                str(self._directory(self.input_path.parent if self.input_path else Path.home()) / "interface_report.json"),
                "JSON files (*.json)",
            )
            if not selected:
                return
            model = self._last_report
            if model is not None and model["config"] != config:
                model = None
            path = Path(selected).expanduser().resolve()
            path.write_text(json.dumps({"config": config, "search": self._result, "model": model},
                                       indent=2, ensure_ascii=False), encoding="utf-8")
            self._remember(path)
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))

    def show_frame(self, metadata: dict[str, Any]) -> None:
        """Follow the main viewer's frame slider without recomputing contacts."""

        if metadata.get("engine") not in {"interface-master", "pymatgen"}:
            self.frame_summary.clear()
            return
        identity = (
            f"XYZ config_type: {metadata['config_type']}\n"
            f"Shift basis: {metadata['shift_basis']} · d1 {metadata['d1']:.1f} · d2 {metadata['d2']:.1f}\n"
            if "config_type" in metadata else ""
        )
        if metadata["engine"] == "pymatgen":
            self.frame_summary.setText(identity.rstrip())
            return
        self.frame_summary.setText(
            identity + f"Current frame {metadata['frame']} · projected coincidence {metadata['registry_fraction']:.1%} "
            f"({metadata['contact_film_atoms']} film contact atoms)\n"
            f"Mean lateral distance {metadata['mean_lateral_distance_A']:.4f} Å · "
            f"minimum cross-interface distance {metadata['minimum_cross_distance_A']:.4f} Å\n"
            f"{metadata['close_contact_count']} film atom(s) below the close-contact distance"
        )

    def choose_structure(self, side: str) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, f"Choose {side} structure",
            str(self._directory(self.input_path.parent if self.input_path else Path.home())), self._file_filter,
        )
        if not selected:
            return
        path = Path(selected).expanduser().resolve()
        self._remember(path)
        if side == "substrate":
            data = copy.deepcopy(self._state)
            data["substrate_structure"] = str(path)
            self._set_config(data)
        else:
            data = copy.deepcopy(self._state)
            data["film_structure"] = str(path)
            self._set_config(data)
        self.changed.emit()

    def load_config(self, *_args: Any) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, "Load interface.json",
            str(self._directory(self.input_path.parent if self.input_path else Path.home())), "JSON files (*.json)",
        )
        if not selected:
            return
        path = Path(selected).expanduser().resolve()
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise ValueError("interface.json must contain a JSON object")
            for key in ("substrate_structure", "film_structure"):
                if data.get(key):
                    source = Path(data[key]).expanduser()
                    data[key] = str((path.parent / source).resolve() if not source.is_absolute() else source.resolve())
            values = self.backend.parameters(data)
            self._invalidate()
            if values.get("substrate_structure"):
                self.input_path = Path(values["substrate_structure"])
                self.substrate_label.setText(str(self.input_path))
            self._set_config(values)
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.config_path = path
        self._remember(path)
        self.changed.emit()

    def save_config(self, *_args: Any) -> None:
        try:
            data = copy.deepcopy(self.config())
        except (TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        if self.input_path is not None:
            data["substrate_structure"] = str(self.input_path)
        selected, _ = QFileDialog.getSaveFileName(
            self, "Save interface.json",
            str(self._directory(self.input_path.parent if self.input_path else Path.home()) / "interface.json"),
            "JSON files (*.json)",
        )
        if selected:
            path = Path(selected).expanduser().resolve()
            try:
                path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            except OSError as exc:
                self.failed.emit(str(exc))
            else:
                self.config_path = path
                self._remember(path)
