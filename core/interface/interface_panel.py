# CCkit 功能说明：
# 异质结建模功能的独立参数卡片，复用主窗口的帧显示和导出。
"""Two-material Interface Builder card for the desktop visualization shell."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Callable

from PySide6.QtCore import QSignalBlocker, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QPlainTextEdit, QPushButton, QVBoxLayout, QWidget,
)


class InterfacePanel(QWidget):
    """Own interface inputs/settings independently of OmniGB's JSON editor."""

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
        self._directory = dialog_directory
        self._remember = remember_directory
        self._file_filter = file_filter
        self._updating = True
        self.numbers: dict[str, Any] = {}
        self.indices: dict[str, Any] = {}
        self.checks: dict[str, Any] = {}
        defaults = backend.default_config()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        files = QGroupBox("Substrate / film")
        form = QFormLayout(files)
        self.substrate_label = QLabel("not selected")
        self.film_label = QLabel("not selected")
        for label in (self.substrate_label, self.film_label):
            label.setWordWrap(True)
        for prefix, label in (("substrate", self.substrate_label), ("film", self.film_label)):
            button = QPushButton(f"Choose {prefix}…")
            button.clicked.connect(lambda _checked=False, side=prefix: self.choose_structure(side))
            form.addRow(label, button)
        layout.addWidget(files)

        contact = QGroupBox("Contact planes and in-plane orientation")
        contact_layout = QVBoxLayout(contact)
        for prefix in ("substrate", "film"):
            row = QHBoxLayout()
            row.addWidget(QLabel(f"{prefix.title()} plane"))
            editor = index_editor("plane")
            editor.set_indices(defaults[f"{prefix}_plane"], "three")
            editor.changed.connect(self.controls_changed)
            editor.invalid.connect(self.failed.emit)
            self.indices[f"{prefix}_plane"] = editor
            row.addWidget(editor)
            contact_layout.addLayout(row)
        self.align_check = QCheckBox("Align specified in-plane reference directions")
        self.align_check.toggled.connect(self.controls_changed)
        contact_layout.addWidget(self.align_check)
        for prefix in ("substrate", "film"):
            row = QHBoxLayout()
            row.addWidget(QLabel(f"{prefix.title()} direction"))
            editor = index_editor("direction")
            editor.set_indices(defaults[f"{prefix}_direction"], "three")
            editor.changed.connect(self.controls_changed)
            editor.invalid.connect(self.failed.emit)
            editor.setEnabled(False)
            self.indices[f"{prefix}_direction"] = editor
            row.addWidget(editor)
            contact_layout.addLayout(row)
        hint = QLabel("Angles are twists about the aligned contact-plane normal.\n"
                      "Each reference direction must lie in its own contact plane.")
        hint.setWordWrap(True)
        contact_layout.addWidget(hint)
        layout.addWidget(contact)

        search = QGroupBox("2D CSL matching")
        search_form = QFormLayout(search)
        self.policy = combo_box()
        self.policy.addItem("Near CSL", "near")
        self.policy.addItem("Exact CSL", "exact")
        self.policy.currentIndexChanged.connect(self.controls_changed)
        search_form.addRow("CSL policy", self.policy)

        def number(form: Any, key: str, label: str, minimum: float, maximum: float,
                   integer: bool = False) -> None:
            spin = integer_spin() if integer else double_spin()
            if integer:
                spin.setRange(int(minimum), int(maximum))
            else:
                spin.setRange(minimum, maximum)
            if not integer:
                spin.setDecimals(8 if key.startswith("tolerance_") else 5)
                spin.setSingleStep(0.01 if "tolerance_" in key else 0.1)
            spin.setKeyboardTracking(False)
            if key.startswith("tolerance_"):
                value = defaults["tolerances"][key.removeprefix("tolerance_")]
            elif key in ("cnid_n1", "cnid_n2"):
                value = defaults["cnid_grid"][0 if key == "cnid_n1" else 1]
            else:
                value = defaults[key]
            spin.setValue(value)
            spin.valueChanged.connect(self.controls_changed)
            self.numbers[key] = spin
            form.addRow(label, spin)

        for key, label in (
            ("angle_start_deg", "Twist start (°)"), ("angle_stop_deg", "Twist stop (°)"),
            ("angle_step_deg", "Twist step (°)"),
        ):
            number(search_form, key, label, 1e-5 if key == "angle_step_deg" else -360, 360)
        for key, label in (
            ("max_sigma_substrate", "Maximum substrate area index"),
            ("max_sigma_film", "Maximum film area index"),
        ):
            number(search_form, key, label, 1, 100000, True)
        for name in ("du", "area", "deformation", "integer"):
            number(search_form, f"tolerance_{name}", f"Tolerance {name}", 1e-8, 1)
        layout.addWidget(search)

        geometry = QGroupBox("Slab geometry and CNID grid")
        geometry_form = QFormLayout(geometry)
        for key, label in (
            ("substrate_normal_repeat", "Substrate normal repeats"),
            ("film_normal_repeat", "Film normal repeats"),
            ("superlattice_lim", "Cell-vector search limit"),
            ("cnid_n1", "CNID u1 subdivisions"), ("cnid_n2", "CNID u2 subdivisions"),
        ):
            number(geometry_form, key, label, 1, 1000, True)
        for key, label in (
            ("dp1", "Substrate termination dp1 (Å)"),
            ("dp2", "Film termination dp2 (Å)"),
            ("normal_distance", "Interface spacing (Å)"), ("vacuum", "External vacuum (Å)"),
        ):
            number(geometry_form, key, label, -1000 if key in {"dp1", "dp2"} else 0, 10000)
        for key, label in (
            ("normal_ortho", "Require normal slab vectors"),
            ("KTI", "Cut / extend slabs when changing termination"),
        ):
            check = QCheckBox(label)
            check.setChecked(defaults[key])
            check.toggled.connect(self.controls_changed)
            self.checks[key] = check
            geometry_form.addRow(check)
        self.summary = QLabel("Select both structures, then build the first admissible match.")
        self.summary.setWordWrap(True)
        geometry_form.addRow(self.summary)
        layout.addWidget(geometry)

        self.auto_check = QCheckBox("Auto refresh")
        self.auto_check.setChecked(True)
        self.auto_check.toggled.connect(lambda _checked: self.changed.emit())
        layout.addWidget(self.auto_check)
        build = QPushButton("Find match and build CNID frames")
        build.clicked.connect(lambda _checked=False: self.build_requested.emit())
        layout.addWidget(build)
        config_group = QGroupBox("interface.json")
        config_layout = QVBoxLayout(config_group)
        self.editor = QPlainTextEdit()
        self.editor.setMinimumHeight(280)
        config_layout.addWidget(self.editor)
        buttons = QHBoxLayout()
        for label, callback in (
            ("Apply JSON", self.apply_json), ("Load JSON…", self.load_config),
            ("Save JSON…", self.save_config),
        ):
            button = QPushButton(label)
            button.clicked.connect(callback)
            buttons.addWidget(button)
        config_layout.addLayout(buttons)
        layout.addWidget(config_group)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(500)
        self._timer.timeout.connect(self.apply_json)
        self.editor.textChanged.connect(self.schedule_json)
        self._set_config(defaults)
        self._updating = False

    def schedule_json(self) -> None:
        if not self._updating:
            self._timer.start()

    def config(self) -> dict[str, Any]:
        """Commit typed numbers and capture the current JSON and film path."""

        for spin in self.numbers.values():
            spin.interpretText()
        for editor in self.indices.values():
            editor.commit()
        data = json.loads(self.editor.toPlainText())
        if not isinstance(data, dict):
            raise ValueError("interface.json must be a JSON object")
        data = self.backend.parameters(data)
        data["mode"] = "build"
        return data

    def controls_changed(self, *_args: Any) -> None:
        if self._updating or not hasattr(self, "editor"):
            return
        try:
            data = json.loads(self.editor.toPlainText())
            if not isinstance(data, dict):
                raise ValueError("interface.json must be a JSON object")
            for key, editor in self.indices.items():
                data[key] = editor.indices()
                data[f"{key}_notation"] = editor.notation()
            data["align_directions"] = self.align_check.isChecked()
            for prefix in ("substrate", "film"):
                self.indices[f"{prefix}_direction"].setEnabled(data["align_directions"])
            data["csl_mode"] = self.policy.currentData()
            data["tolerances"] = {
                name: self.numbers[f"tolerance_{name}"].value()
                for name in ("du", "area", "deformation", "integer")
            }
            data["cnid_grid"] = [self.numbers["cnid_n1"].value(), self.numbers["cnid_n2"].value()]
            for key, spin in self.numbers.items():
                if not key.startswith("tolerance_") and key not in {"cnid_n1", "cnid_n2"}:
                    data[key] = spin.value()
            data.update({key: check.isChecked() for key, check in self.checks.items()})
            self._updating = True
            self.editor.setPlainText(json.dumps(data, indent=2, ensure_ascii=False))
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            self.failed.emit(str(exc))
        else:
            self.changed.emit()
        finally:
            self._updating = False

    def _set_config(self, config: dict[str, Any]) -> None:
        values = self.backend.parameters(config)
        self._updating = True
        try:
            self.editor.setPlainText(json.dumps(values, indent=2, ensure_ascii=False))
            for key, editor in self.indices.items():
                editor.set_indices(values[key], values[f"{key}_notation"])
                if key.endswith("_direction"):
                    editor.setEnabled(values["align_directions"])
            self.align_check.setChecked(values["align_directions"])
            for key, spin in self.numbers.items():
                if key.startswith("tolerance_"):
                    value = values["tolerances"][key.removeprefix("tolerance_")]
                elif key in {"cnid_n1", "cnid_n2"}:
                    value = values["cnid_grid"][0 if key == "cnid_n1" else 1]
                else:
                    value = values[key]
                blocker = QSignalBlocker(spin)
                spin.setValue(value)
                del blocker
            for key, check in self.checks.items():
                check.setChecked(values[key])
            self.policy.setCurrentIndex(self.policy.findData(values["csl_mode"]))
            self.film_label.setText(str(values.get("film_structure") or "not selected"))
        finally:
            self._updating = False

    def apply_json(self, *_args: Any) -> None:
        try:
            self._set_config(self.config())
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            self.failed.emit(str(exc))
            return
        self.changed.emit()

    def choose_structure(self, side: str) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, f"Choose {side} structure",
            str(self._directory(self.input_path.parent if self.input_path else Path.home())),
            self._file_filter,
        )
        if not selected:
            return
        path = Path(selected).expanduser().resolve()
        self._remember(path)
        if side == "substrate":
            self.input_path = path
            self.substrate_label.setText(str(path))
            self.changed.emit()
        else:
            try:
                config = self.config()
                config["film_structure"] = str(path)
                self._set_config(config)
            except (ValueError, TypeError) as exc:
                self.failed.emit(str(exc))
                return
            self.changed.emit()

    def load_config(self, *_args: Any) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, "Load interface.json",
            str(self._directory(self.input_path.parent if self.input_path else Path.home())),
            "JSON files (*.json)",
        )
        if not selected:
            return
        path = Path(selected).expanduser().resolve()
        try:
            config = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(config, dict):
                raise ValueError("interface.json must be a JSON object")
            for key in ("substrate_structure", "film_structure"):
                if config.get(key):
                    source = Path(config[key]).expanduser()
                    if not source.is_absolute():
                        source = path.parent / source
                    config[key] = str(source.resolve())
            self._set_config(config)
            if config.get("substrate_structure"):
                self.input_path = Path(config["substrate_structure"])
                self.substrate_label.setText(str(self.input_path))
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.config_path = path
        self._remember(path)
        self.changed.emit()

    def save_config(self, *_args: Any) -> None:
        try:
            config = copy.deepcopy(self.config())
        except (TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        if self.input_path is not None:
            config["substrate_structure"] = str(self.input_path)
        directory = self.input_path.parent if self.input_path else Path.home()
        selected, _ = QFileDialog.getSaveFileName(
            self, "Save interface.json", str(directory / "interface.json"), "JSON files (*.json)",
        )
        if selected:
            path = Path(selected).expanduser().resolve()
            try:
                path.write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            except OSError as exc:
                self.failed.emit(str(exc))
            else:
                self.config_path = path
                self._remember(path)
