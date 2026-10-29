"""Independent crystal orientation, slab geometry and termination controls."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QPlainTextEdit, QPushButton, QSizePolicy, QVBoxLayout, QWidget,
)

from viewer_widgets import WrappingLabel as QLabel


class SurfacePanel(QWidget):
    """Own Surface Builder input and settings across feature switches."""

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
        self._orientation_key: str | None = None
        self._choices: list[dict[str, Any]] = []
        self.indices: dict[str, Any] = {}
        self.numbers: dict[str, Any] = {}
        self.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        def group(title: str) -> QFormLayout:
            box = QGroupBox(title)
            box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
            form = QFormLayout(box)
            form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
            form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
            layout.addWidget(box)
            return form

        def number(form: QFormLayout, key: str, label: str, lower: float, upper: float,
                   *, integer: bool = False, tip: str = "") -> None:
            spin = integer_spin() if integer else double_spin()
            spin.setRange(int(lower), int(upper)) if integer else spin.setRange(lower, upper)
            if not integer:
                spin.setDecimals(6)
                spin.setSingleStep(0.01 if key in {"cut_phase", "layer_tolerance"} else 0.5)
            spin.setKeyboardTracking(False)
            spin.setToolTip(tip)
            spin.valueChanged.connect(self._edited)
            self.numbers[key] = spin
            form.addRow(label, spin)

        inputs = group("1 · Bulk structure and surface orientation")
        self.input_label = QLabel("not selected")
        self.input_label.setWordWrap(True)
        self.input_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        choose = QPushButton("Choose bulk…")
        choose.clicked.connect(self.choose_structure)
        inputs.addRow(choose, self.input_label)
        for key, kind, label in (
            ("plane", "plane", "Surface plane"),
            ("direction_a", "direction", "In-plane direction a"),
            ("direction_b", "direction", "In-plane direction b"),
        ):
            editor = index_editor(kind)
            editor.changed.connect(self._edited)
            editor.invalid.connect(self.failed.emit)
            self.indices[key] = editor
            inputs.addRow(label, editor)
        note = QLabel(
            "Indices refer to the imported bulk cell. Both directions must lie in the plane "
            "(hU + kV + lW = 0) and form a right-handed pair with its normal. "
            "Miller–Bravais indices require a standard hexagonal cell."
        )
        note.setWordWrap(True)
        inputs.addRow(note)
        number(inputs, "repeat_a", "Repeats along a", 1, 100, integer=True)
        number(inputs, "repeat_b", "Repeats along b", 1, 100, integer=True)

        geometry = group("2 · Thickness and vacuum")
        number(geometry, "thickness", "Cut-window thickness (Å)", 0.001, 10000,
               tip="Normal width of the bulk interval to retain. The outermost-atom thickness is discrete and may be smaller.")
        number(geometry, "vacuum", "Vacuum gap (Å)", 0.001, 10000,
               tip="Normal distance between the outermost atomic planes of neighboring periodic slabs; split equally above and below this slab.")
        number(geometry, "max_atoms", "Maximum atoms", 1, 200000, integer=True)
        geometry_note = QLabel(
            "The unstrained slab is centered in vacuum. Output a/b follow the two chosen "
            "directions; c is perpendicular to the surface."
        )
        geometry_note.setWordWrap(True)
        geometry.addRow(geometry_note)

        termination = group("3 · Surface termination")
        self.termination_mode = combo_box()
        self.termination_mode.addItem("Select an interlayer cut", "layer")
        self.termination_mode.addItem("Custom cut position", "custom")
        self.termination_mode.currentIndexChanged.connect(self._edited)
        termination.addRow("Cut mode", self.termination_mode)
        self.layer_panel = QWidget()
        layer_form = QFormLayout(self.layer_panel)
        layer_form.setContentsMargins(0, 0, 0, 0)
        layer_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self.cut_combo = combo_box()
        self.cut_combo.setMinimumContentsLength(12)
        self.cut_combo.setSizeAdjustPolicy(self.cut_combo.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.cut_combo.currentIndexChanged.connect(self._edited)
        layer_form.addRow("Interlayer cut", self.cut_combo)
        termination.addRow(self.layer_panel)
        self.custom_panel = QWidget()
        custom_form = QFormLayout(self.custom_panel)
        custom_form.setContentsMargins(0, 0, 0, 0)
        number(custom_form, "cut_phase", "Bottom cut phase (fraction)", 0, 1,
               tip="Lower boundary along the positive surface normal, in units of the reduced (hkl) plane spacing. 1 and 0 are the same phase.")
        termination.addRow(self.custom_panel)
        number(termination, "layer_tolerance", "Layer tolerance (Å)", 0.000001, 5,
               tip="Gaps smaller than this value are treated as belonging to one atomic layer. Also used to identify the exposed elements.")
        self.termination_hint = QLabel()
        self.termination_hint.setWordWrap(True)
        termination.addRow(self.termination_hint)

        self.auto_check = QCheckBox("Auto rebuild after parameter changes")
        self.auto_check.toggled.connect(self._edited)
        layout.addWidget(self.auto_check)
        build = QPushButton("Generate surface")
        build.clicked.connect(self._generate)
        layout.addWidget(build)
        self.summary = QLabel("Choose a bulk structure and orientation, then generate the surface.")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.summary)
        self.vectors = QPlainTextEdit()
        self.vectors.setReadOnly(True)
        self.vectors.setMinimumHeight(85)
        self.vectors.setMaximumHeight(130)
        layout.addWidget(self.vectors)
        buttons = QHBoxLayout()
        for label, callback in (("Load JSON…", self.load_config), ("Save JSON…", self.save_config)):
            button = QPushButton(label)
            button.clicked.connect(callback)
            buttons.addWidget(button)
        layout.addLayout(buttons)
        self._set_config(backend.default_config())
        self._updating = False

    def _values(self) -> dict[str, Any]:
        values = {"bulk_structure": str(self.input_path) if self.input_path else ""}
        for key, editor in self.indices.items():
            values[key], values[f"{key}_notation"] = editor.indices(), editor.notation()
        values.update({key: spin.value() for key, spin in self.numbers.items()})
        values["termination_mode"] = self.termination_mode.currentData()
        values["termination_index"] = self.cut_combo.currentData() or 0
        return self.backend.parameters(values)

    def _key(self, values: dict[str, Any]) -> str:
        return json.dumps({key: values[key] for key in (
            "bulk_structure", "plane", "plane_notation", "direction_a", "direction_a_notation",
            "direction_b", "direction_b_notation", "repeat_a", "repeat_b", "layer_tolerance",
        )}, sort_keys=True)

    def _show_modes(self) -> None:
        custom = self.termination_mode.currentData() == "custom"
        self.layer_panel.setVisible(not custom)
        self.custom_panel.setVisible(custom)
        self.termination_hint.setText(
            "Move the lower cutting plane through one normal repeat; the upper cut follows "
            "at the selected thickness. This changes which atomic layers are exposed."
            if custom else
            "Generate once to list all interlayer cuts. Choose another cut to change the "
            "exposed layers. Top and bottom terminations are coupled to the thickness."
        )

    def _edited(self, *_args: Any) -> None:
        if self._updating:
            return
        self._show_modes()
        try:
            values = self._values()
            key = self._key(values)
        except (TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            self.changed.emit()
            return
        if key != self._orientation_key:
            self._updating = True
            self._choices = []
            self.cut_combo.clear()
            self.cut_combo.addItem("Cut 1 · generate to list layers", 0)
            self._orientation_key = key
            self._updating = False
        self.summary.setText("Parameters changed; generate the surface to update the model.")
        self.vectors.clear()
        self.changed.emit()

    def config(self) -> dict[str, Any]:
        """Commit active editors and capture a normalized build request."""

        for editor in self.indices.values():
            editor.commit()
        for spin in self.numbers.values():
            spin.interpretText()
        return self._values()

    def _set_config(self, config: dict[str, Any]) -> None:
        values = self.backend.parameters(config)
        previous = self._updating
        self._updating = True
        try:
            self.input_path = Path(values["bulk_structure"]).expanduser().resolve() if values["bulk_structure"] else None
            self.input_label.setText(str(self.input_path) if self.input_path else "not selected")
            for key, editor in self.indices.items():
                editor.set_indices(values[key], values[f"{key}_notation"])
            for key, spin in self.numbers.items():
                spin.setValue(values[key])
            self.termination_mode.setCurrentIndex(self.termination_mode.findData(values["termination_mode"]))
            self._choices = []
            self.cut_combo.clear()
            index = values["termination_index"]
            self.cut_combo.addItem(f"Cut {index + 1} · generate to list layers", index)
            self._orientation_key = self._key(self._values())
            self._show_modes()
            self.summary.setText("Generate the surface for the current parameters.")
            self.vectors.clear()
        finally:
            self._updating = previous

    def _generate(self) -> None:
        try:
            self.config()
            if self.input_path is None or not self.input_path.is_file():
                raise ValueError("Choose a bulk structure before generating a surface")
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.build_requested.emit()

    def choose_structure(self) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, "Choose bulk structure", str(self._directory(self.input_path.parent if self.input_path else Path.cwd())),
            self._file_filter,
        )
        if selected:
            self.input_path = Path(selected).resolve()
            self._remember(self.input_path)
            self.input_label.setText(str(self.input_path))
            self._edited()

    def show_report(self, report: dict[str, Any]) -> None:
        """List available cuts without scheduling a second build from result updates."""

        self._updating = True
        try:
            selected = self.cut_combo.currentData()
            self._choices = report["cut_choices"]
            self.cut_combo.clear()
            for index, choice in enumerate(self._choices):
                label = ", ".join(choice["bottom_species"])
                text = f"Cut {index + 1} · phase {choice['phase']:.6f} · bottom {label}"
                self.cut_combo.addItem(text, index)
                self.cut_combo.setItemData(index, text, Qt.ItemDataRole.ToolTipRole)
            if selected is not None and self.cut_combo.findData(selected) < 0:
                # A custom cut does not use the saved interlayer selection.
                self.cut_combo.addItem(f"Cut {selected + 1} · unavailable for current layers", selected)
            self.cut_combo.setCurrentIndex(max(0, self.cut_combo.findData(selected)))
        finally:
            self._updating = False
        self.summary.setText(
            f"{report['atom_count']} atoms · surface area {report['area_A2']:.3f} Å²\n"
            f"Cut window {report['cut_window_A']:.3f} Å · atomic thickness {report['atomic_thickness_A']:.3f} Å\n"
            f"Vacuum gap {report['vacuum_A']:.3f} Å · cut phase {report['cut_phase']:.6f}\n"
            f"Top layer: {', '.join(report['top_species'])} · bottom layer: {', '.join(report['bottom_species'])}\n"
            f"{len(self._choices)} interlayer cuts · normal repeat {report['normal_period_A']:.4f} Å"
        )
        lines = ["In-plane vectors a / b (Å):"]
        lines.extend("  " + "  ".join(f"{item:.6f}" for item in row) for row in report["in_plane_vectors_A"])
        lines.append("Oriented-cell matrix (input bulk basis):")
        lines.extend("  " + str(row) for row in report["oriented_matrix"])
        self.vectors.setPlainText("\n".join(lines))

    def load_config(self) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, "Load surface JSON", str(self._directory(self.input_path.parent if self.input_path else Path.cwd())),
            "JSON files (*.json)",
        )
        if not selected:
            return
        try:
            path = Path(selected).resolve()
            values = self.backend.parameters(json.loads(path.read_text(encoding="utf-8")))
            if values["bulk_structure"]:
                source = Path(values["bulk_structure"]).expanduser()
                values["bulk_structure"] = str(source if source.is_absolute() else path.parent / source)
            self._set_config(values)
            self.config_path = path
            self._remember(path)
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.changed.emit()

    def save_config(self) -> None:
        try:
            values = self.config()
            selected, _ = QFileDialog.getSaveFileName(
                self, "Save surface JSON", str(self.config_path or self._directory(
                    self.input_path.parent if self.input_path else Path.cwd(),
                ) / "surface.json"), "JSON files (*.json)",
            )
            if selected:
                path = Path(selected).resolve()
                path.write_text(json.dumps(values, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
                self.config_path = path
                self._remember(path)
        except (OSError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
