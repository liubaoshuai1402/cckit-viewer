# CCkit 功能说明：Wyckoff 晶体建模的独立参数卡片。
"""Space-group, lattice and editable Wyckoff-site controls for the viewer."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

from PySide6.QtCore import QSignalBlocker, Qt, Signal
from PySide6.QtWidgets import (
    QFileDialog, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QSizePolicy, QVBoxLayout, QWidget,
)


class WyckoffPanel(QWidget):
    """Own a manually generated crystal independently of file-based modelers."""

    changed = Signal()
    build_requested = Signal()
    failed = Signal(str)

    def __init__(
        self, backend: Any, integer_spin: Any, double_spin: Any, combo_box: Any,
        dialog_directory: Callable[[Path], Path], remember_directory: Callable[[Path], None],
    ) -> None:
        super().__init__()
        self.backend = backend
        self._combo_box = combo_box
        self._directory = dialog_directory
        self._remember = remember_directory
        self.config_path: Path | None = None
        self.rows: list[dict[str, Any]] = []
        self._updating = True
        self.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        card = QGroupBox("wyckoff model")
        card.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(10, 14, 10, 10)
        card_layout.setSpacing(8)
        card_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        group_row = QHBoxLayout()
        group_row.setSpacing(8)
        group_row.addWidget(QLabel("Space group"))
        self.group_number = integer_spin()
        self.group_number.setRange(1, 230)
        self.group_number.setValue(194)
        self.group_number.setMaximumWidth(100)
        self.group_number.setKeyboardTracking(False)
        self.group_number.valueChanged.connect(self._group_changed)
        group_row.addWidget(self.group_number)
        self.symbol = QLabel()
        self.symbol.setStyleSheet("font-weight: 600; color: #334155;")
        group_row.addWidget(self.symbol, 1)
        card_layout.addLayout(group_row)
        self.setting_label = QLabel()
        self.setting_label.setWordWrap(True)
        self.setting_label.setStyleSheet("color: #64748b;")
        card_layout.addWidget(self.setting_label)

        self.lattice_group = QGroupBox("Lattice constants (Å); angles (°)")
        self.lattice_group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        self.lattice_layout = QGridLayout(self.lattice_group)
        self.lattice_layout.setContentsMargins(10, 14, 10, 10)
        self.lattice_layout.setHorizontalSpacing(10)
        self.lattice_layout.setVerticalSpacing(6)
        self.lattice_layout.setColumnStretch(1, 1)
        self.lattice_controls: dict[str, Any] = {}
        self.lattice_labels: dict[str, Any] = {}
        for index, key in enumerate(("a", "b", "c", "alpha", "beta", "gamma")):
            label = QLabel(key + " =", self.lattice_group)
            spin = double_spin(self.lattice_group)
            spin.setDecimals(6)
            spin.setRange(0.000001, 100000.0 if index < 3 else 179.999999)
            spin.setValue(3.231 if index < 3 else 90.0)
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(self._edited)
            self.lattice_labels[key] = label
            self.lattice_controls[key] = spin
        self.lattice_controls["c"].setValue(5.147)
        card_layout.addWidget(self.lattice_group)

        sites_group = QGroupBox("Wyckoff positions")
        sites_group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        self.sites_layout = QVBoxLayout(sites_group)
        self.sites_layout.setContentsMargins(10, 14, 10, 10)
        self.sites_layout.setSpacing(8)
        self.sites_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        hint = QLabel("Choose equivalent fractional coordinates, or select Manual input.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #64748b;")
        self.sites_layout.addWidget(hint)
        card_layout.addWidget(sites_group)
        self.add_button = QPushButton("Add Wyckoff position")
        self.add_button.clicked.connect(lambda: self.add_site())
        card_layout.addWidget(self.add_button)
        self.generate_button = QPushButton("Generate crystal")
        self.generate_button.setStyleSheet(
            "QPushButton { background: #2563eb; color: white; border: none; "
            "border-radius: 4px; padding: 7px 12px; font-weight: 600; }"
            "QPushButton:hover { background: #1d4ed8; }"
            "QPushButton:pressed { background: #1e40af; }"
            "QPushButton:disabled { background: #94a3b8; }"
        )
        self.generate_button.clicked.connect(self._generate)
        card_layout.addWidget(self.generate_button)
        self.summary = QLabel("Fill in the parameters, then generate the crystal.")
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet("color: #475569;")
        card_layout.addWidget(self.summary)
        config_buttons = QHBoxLayout()
        for title, callback in (("Load JSON…", self.load_config), ("Save JSON…", self.save_config)):
            button = QPushButton(title)
            button.clicked.connect(callback)
            config_buttons.addWidget(button)
        card_layout.addLayout(config_buttons)
        layout.addWidget(card)
        self._group_changed()
        if self.generate_button.isEnabled():
            self.add_site(backend.default_config()["sites"][0])
        self._updating = False

    @property
    def input_path(self) -> Path:
        """Supply a stable request identity; this feature needs no input file."""

        return self.config_path or Path.home() / "wyckoff_model.cif"

    def _group_changed(self, *_args: Any) -> None:
        number = self.group_number.value()
        try:
            group = self.backend.space_group(number)
        except (ImportError, RuntimeError, ValueError) as exc:
            self.symbol.setText("PyXtal unavailable")
            self.setting_label.setText(str(exc))
            self.add_button.setEnabled(False)
            self.generate_button.setEnabled(False)
            self.failed.emit(str(exc))
            return
        self.generate_button.setEnabled(True)
        self.add_button.setEnabled(True)
        self.symbol.setText(group.symbol)
        system = group.lattice_type.title()
        setting = "Conventional cell · spglib setting"
        if 3 <= number <= 15:
            setting += " · unique axis b"
        elif 143 <= number <= 194:
            setting += " · hexagonal axes, a=b, γ=120°"
        self.setting_label.setText(f"{system} · {setting}")
        fields = self.backend.lattice_fields(number)
        for key, spin in self.lattice_controls.items():
            self.lattice_layout.removeWidget(spin)
            self.lattice_layout.removeWidget(self.lattice_labels[key])
            spin.setVisible(key in fields)
            self.lattice_labels[key].setVisible(key in fields)
        for index, key in enumerate(fields):
            self.lattice_layout.addWidget(self.lattice_labels[key], index, 0)
            self.lattice_layout.addWidget(self.lattice_controls[key], index, 1)
        for row in self.rows:
            self._fill_positions(row)
        self._edited()

    def _fill_positions(self, row: dict[str, Any], selected: str | None = None) -> None:
        selector = row["selector"]
        blocker = QSignalBlocker(selector)
        selector.clear()
        for position in self.backend.space_group(self.group_number.value()):
            selector.addItem(position.get_label(), position.get_label())
        if selected is not None:
            selector.setCurrentIndex(max(0, selector.findData(selected)))
        del blocker
        self._position_changed(row)

    def add_site(self, site: dict[str, Any] | None = None) -> None:
        """Append an independently editable Wyckoff orbit."""

        widget = QFrame()
        widget.setObjectName("wyckoffSiteCard")
        widget.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        widget.setStyleSheet(
            "QFrame#wyckoffSiteCard { background: #ffffff; border: 1px solid #cbd5e1; "
            "border-radius: 6px; }"
        )
        site_layout = QVBoxLayout(widget)
        site_layout.setContentsMargins(10, 8, 10, 8)
        site_layout.setSpacing(6)
        header = QHBoxLayout()
        title = QLabel()
        title.setStyleSheet("color: #334155; font-weight: 600;")
        header.addWidget(title)
        header.addStretch(1)
        remove = QPushButton("Remove")
        remove.setMaximumWidth(90)
        header.addWidget(remove)
        site_layout.addLayout(header)
        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(6)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)
        element = QLineEdit((site or {}).get("element", ""))
        element.setMinimumWidth(0)
        element.setPlaceholderText("e.g. Zr")
        selector = self._combo_box()
        grid.addWidget(QLabel("Element"), 0, 0)
        grid.addWidget(element, 0, 1)
        grid.addWidget(QLabel("Position"), 0, 2)
        grid.addWidget(selector, 0, 3)
        site_layout.addLayout(grid)
        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Coordinates"))
        coordinate_mode = self._combo_box()
        coordinate_mode.addItem("Select equivalent coordinates", "equivalent")
        coordinate_mode.addItem("Manual input", "manual")
        mode_row.addWidget(coordinate_mode, 1)
        site_layout.addLayout(mode_row)

        automatic_widget = QWidget()
        automatic_layout = QVBoxLayout(automatic_widget)
        automatic_layout.setContentsMargins(0, 0, 0, 0)
        automatic_layout.setSpacing(4)
        equivalent_selector = self._combo_box()
        equivalent_selector.setMinimumWidth(0)
        equivalent_selector.setMinimumContentsLength(12)
        equivalent_selector.setSizeAdjustPolicy(
            equivalent_selector.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon,
        )
        automatic_layout.addWidget(equivalent_selector)
        parameter_group = QGroupBox("Free parameters (fractional)")
        parameter_group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        parameter_layout = QGridLayout(parameter_group)
        parameter_layout.setContentsMargins(8, 14, 8, 8)
        parameter_layout.setHorizontalSpacing(8)
        parameter_layout.setVerticalSpacing(6)
        parameter_layout.setColumnStretch(1, 1)
        parameter_labels = {}
        parameter_editors = {}
        for name in ("x", "y", "z"):
            parameter_labels[name] = QLabel(name + " =", parameter_group)
            editor = QLineEdit(parameter_group)
            editor.setMinimumWidth(0)
            editor.setPlaceholderText("Decimal or fraction, e.g. 1/3")
            parameter_editors[name] = editor
        parameter_error = QLabel(parameter_group)
        parameter_error.setWordWrap(True)
        parameter_error.setStyleSheet("color: #b42318;")
        parameter_error.hide()
        parameter_layout.addWidget(parameter_error, 3, 0, 1, 2)
        automatic_layout.addWidget(parameter_group)
        choice_note = QLabel()
        choice_note.setWordWrap(True)
        choice_note.setStyleSheet("color: #64748b;")
        automatic_layout.addWidget(choice_note)
        site_layout.addWidget(automatic_widget)

        manual_widget = QWidget()
        manual_layout = QVBoxLayout(manual_widget)
        manual_layout.setContentsMargins(0, 0, 0, 0)
        manual_layout.setSpacing(4)
        coordinate_row = QHBoxLayout()
        coordinate_row.setSpacing(8)
        coordinates = []
        for name in ("x", "y", "z"):
            editor = QLineEdit()
            editor.setMinimumWidth(0)
            editor.setPlaceholderText(name)
            coordinate_row.addWidget(QLabel(name + " ="))
            coordinate_row.addWidget(editor, 1)
            coordinates.append(editor)
        manual_layout.addLayout(coordinate_row)
        manual_hint = QLabel(
            "Fractional x, y, z · decimals or fractions (e.g. 1/3). "
            "Switch to Select equivalent coordinates to list this orbit."
        )
        manual_hint.setWordWrap(True)
        manual_hint.setStyleSheet("color: #64748b;")
        manual_layout.addWidget(manual_hint)
        site_layout.addWidget(manual_widget)
        constraint = QLabel()
        constraint.setWordWrap(True)
        constraint.setStyleSheet("color: #64748b;")
        site_layout.addWidget(constraint)
        row = {
            "widget": widget, "element": element, "selector": selector,
            "coordinates": coordinates, "constraint": constraint, "title": title,
            "coordinate_mode": coordinate_mode, "equivalent_selector": equivalent_selector,
            "automatic_widget": automatic_widget, "manual_widget": manual_widget,
            "choice_note": choice_note, "updating": False, "example": True,
            "parameter_group": parameter_group, "parameter_layout": parameter_layout,
            "parameter_labels": parameter_labels, "parameter_editors": parameter_editors,
            "parameter_names": [], "parameters_dirty": False, "parameter_error": parameter_error,
        }
        self.rows.append(row)
        self._update_site_titles()
        self.sites_layout.addWidget(widget)
        selector.currentIndexChanged.connect(lambda: self._position_changed(row))
        coordinate_mode.currentIndexChanged.connect(lambda: self._coordinate_mode_changed(row))
        equivalent_selector.currentIndexChanged.connect(lambda: self._equivalent_changed(row))
        for editor in coordinates:
            editor.textChanged.connect(lambda _text, active=row: self._manual_coordinates_changed(active))
        for editor in parameter_editors.values():
            editor.textChanged.connect(lambda _text, active=row: self._parameter_text_changed(active))
            editor.editingFinished.connect(lambda active=row: self._parameters_changed(active))
        element.textChanged.connect(self._edited)
        remove.clicked.connect(lambda: self._remove_site(row))
        self._fill_positions(row, (site or {}).get("wyckoff"))
        if site:
            self._restore_site_coordinates(row, site)
        self._edited()

    def _position_changed(self, row: dict[str, Any]) -> None:
        label = row["selector"].currentData()
        options = self.backend.position_options(self.group_number.value(), label)
        row["example"] = bool(options["free_parameters"])
        self._set_equivalent_options(row, options)
        row["updating"] = True
        try:
            row["coordinate_mode"].setCurrentIndex(0)
        finally:
            row["updating"] = False
        self._show_coordinate_mode(row)
        self._edited()

    @staticmethod
    def _format_coordinate(value: float) -> str:
        """Display exact crystallographic fractions and retain generic decimals."""

        fraction = Fraction(float(value)).limit_denominator(96)
        if abs(float(fraction) - value) < 1e-9:
            return str(fraction)
        return f"{value:.12g}"

    def _set_equivalent_options(
        self, row: dict[str, Any], options: dict[str, Any], selected: list[Any] | None = None,
    ) -> None:
        """Fill every member of one orbit and preserve a supplied representative."""

        row["updating"] = True
        try:
            selector = row["equivalent_selector"]
            selector.clear()
            preferred = (
                [self.backend.coordinate_value(value) % 1.0 for value in selected]
                if selected is not None else None
            )
            selected_index = 0
            for index, equivalent in enumerate(options["equivalents"]):
                coords = equivalent["coords"]
                numeric = ", ".join(self._format_coordinate(value) for value in coords)
                text = f"{index + 1}. ({numeric})"
                if options["free_parameters"]:
                    text = f"{index + 1}. ({equivalent['expression']}) → ({numeric})"
                selector.addItem(text, coords)
                selector.setItemData(index, text, Qt.ItemDataRole.ToolTipRole)
                if preferred is not None and all(
                    abs((first - second + 0.5) % 1.0 - 0.5) < 1e-5
                    for first, second in zip(preferred, coords)
                ):
                    selected_index = index
            selector.setCurrentIndex(selected_index)
            row["constraint"].setText(
                f"{options['label']}: ({options['constraint']}) · "
                f"{options['multiplicity']} equivalent atoms"
            )
            free_count = options["free_parameters"]
            row["parameter_names"] = options["parameter_names"]
            for name, editor in row["parameter_editors"].items():
                label = row["parameter_labels"][name]
                row["parameter_layout"].removeWidget(label)
                row["parameter_layout"].removeWidget(editor)
                visible = name in row["parameter_names"]
                label.setVisible(visible)
                editor.setVisible(visible)
                if visible:
                    editor.setText(self._format_coordinate(options["parameter_values"][name]))
            for index, name in enumerate(row["parameter_names"]):
                row["parameter_layout"].addWidget(row["parameter_labels"][name], index, 0)
                row["parameter_layout"].addWidget(row["parameter_editors"][name], index, 1)
            row["parameter_group"].setVisible(bool(row["parameter_names"]))
            row["parameters_dirty"] = False
            row["parameter_error"].clear()
            row["parameter_error"].hide()
            if not free_count:
                note = "Fixed position. Every listed representative generates the same complete orbit."
            elif row["example"]:
                note = (
                    f"{free_count} free parameter(s): this list shows one example orbit. "
                    "Set the parameters above, or use Manual input; the parameters have infinitely many values."
                )
            else:
                note = (
                    f"{free_count} free parameter(s). Listed coordinates form the complete orbit "
                    "for the supplied parameters. Parameter values can vary continuously."
                )
            row["choice_note"].setText(note)
        finally:
            row["updating"] = False
        self._equivalent_changed(row)

    def _equivalent_changed(self, row: dict[str, Any]) -> None:
        if row["updating"]:
            return
        selector = row["equivalent_selector"]
        coords = selector.currentData()
        if coords is None:
            return
        row["updating"] = True
        try:
            for editor, value in zip(row["coordinates"], coords):
                editor.setText(self._format_coordinate(value))
            selector.setToolTip(selector.currentText())
        finally:
            row["updating"] = False
        self._edited()

    def _parameter_text_changed(self, row: dict[str, Any]) -> None:
        if row["updating"]:
            return
        row["parameters_dirty"] = True
        row["example"] = False
        self._edited()

    def _parameters_changed(self, row: dict[str, Any]) -> None:
        """Recalculate the orbit on commit, retaining the selected expression."""

        if row["updating"] or not row["parameters_dirty"]:
            return
        if row["coordinate_mode"].currentData() != "equivalent":
            return
        supplied = {name: row["parameter_editors"][name].text() for name in row["parameter_names"]}
        try:
            options = self.backend.position_options(
                self.group_number.value(), row["selector"].currentData(), free_values=supplied,
            )
        except (RuntimeError, TypeError, ValueError) as exc:
            row["parameter_error"].setText(str(exc))
            row["parameter_error"].show()
            self._edited()
            return
        index = max(0, row["equivalent_selector"].currentIndex())
        self._set_equivalent_options(row, options, options["equivalents"][index]["coords"])
        self._edited()

    def _manual_coordinates_changed(self, row: dict[str, Any]) -> None:
        if not row["updating"]:
            row["example"] = False
            self._edited()

    def _show_coordinate_mode(self, row: dict[str, Any]) -> None:
        manual = row["coordinate_mode"].currentData() == "manual"
        row["manual_widget"].setVisible(manual)
        row["automatic_widget"].setVisible(not manual)
        self._edited()

    def _coordinate_mode_changed(self, row: dict[str, Any]) -> None:
        if row["updating"]:
            return
        if row["coordinate_mode"].currentData() == "equivalent":
            supplied = [editor.text() for editor in row["coordinates"]]
            try:
                options = self.backend.position_options(
                    self.group_number.value(), row["selector"].currentData(), supplied,
                )
            except (RuntimeError, TypeError, ValueError) as exc:
                blocker = QSignalBlocker(row["coordinate_mode"])
                row["coordinate_mode"].setCurrentIndex(1)
                del blocker
                self.failed.emit(f"Could not list equivalent coordinates: {exc}")
            else:
                self._set_equivalent_options(row, options, supplied)
        self._show_coordinate_mode(row)

    def _restore_site_coordinates(self, row: dict[str, Any], site: dict[str, Any]) -> None:
        """Restore mode and coordinates, treating older JSON files as manual."""

        mode = site.get("coordinate_mode", "manual")
        supplied = site["coords"]
        row["example"] = False
        if mode == "equivalent":
            try:
                options = self.backend.position_options(
                    self.group_number.value(), row["selector"].currentData(), supplied,
                )
            except (RuntimeError, TypeError, ValueError) as exc:
                mode = "manual"
                self.failed.emit(f"Loaded coordinates require correction: {exc}")
            else:
                self._set_equivalent_options(row, options, supplied)
        row["updating"] = True
        try:
            if mode == "manual":
                for editor, value in zip(row["coordinates"], supplied):
                    editor.setText(str(value))
            row["coordinate_mode"].setCurrentIndex(row["coordinate_mode"].findData(mode))
        finally:
            row["updating"] = False
        self._show_coordinate_mode(row)

    def _site_config(self, row: dict[str, Any]) -> dict[str, Any]:
        mode = row["coordinate_mode"].currentData()
        if mode == "equivalent":
            self._parameters_changed(row)
            if row["parameters_dirty"]:
                raise ValueError(f"{row['title'].text()}: {row['parameter_error'].text()}")
        coords = (
            [editor.text() for editor in row["coordinates"]]
            if mode == "manual" else row["equivalent_selector"].currentData()
        )
        return {
            "element": row["element"].text(), "wyckoff": row["selector"].currentData(),
            "coords": coords, "coordinate_mode": mode,
        }

    def _remove_site(self, row: dict[str, Any]) -> None:
        self.rows.remove(row)
        self.sites_layout.removeWidget(row["widget"])
        row["widget"].deleteLater()
        self._update_site_titles()
        self._edited()

    def _update_site_titles(self) -> None:
        """Keep the small orbit cards numbered after insertion or removal."""

        for index, row in enumerate(self.rows, 1):
            row["title"].setText(f"Site {index}")

    def _edited(self, *_args: Any) -> None:
        self.updateGeometry()
        if not self._updating:
            self.summary.setText("Parameters changed. Click Generate crystal to update the preview.")
            self.changed.emit()

    def config(self) -> dict[str, Any]:
        """Commit controls and return validated fractional-coordinate inputs."""

        self.group_number.interpretText()
        fields = self.backend.lattice_fields(self.group_number.value())
        for key in fields:
            self.lattice_controls[key].interpretText()
        return self.backend.parameters({
            "space_group": self.group_number.value(),
            "lattice": {key: self.lattice_controls[key].value() for key in fields},
            "sites": [self._site_config(row) for row in self.rows],
        })

    def _generate(self) -> None:
        try:
            self.config()
        except (RuntimeError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self.build_requested.emit()

    def load_config(self, *_args: Any) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, "Load wyckoff.json", str(self._directory(self.input_path.parent)), "JSON files (*.json)",
        )
        if not selected:
            return
        path = Path(selected).expanduser().resolve()
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise ValueError("wyckoff.json must contain a JSON object")
            values = self.backend.parameters(data)
            for key, value in values["lattice"].items():
                spin = self.lattice_controls[key]
                if not spin.minimum() <= value <= spin.maximum():
                    raise ValueError(f"Lattice parameter {key} exceeds the GUI input range")
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        self._updating = True
        try:
            for row in self.rows[:]:
                self._remove_site(row)
            self.group_number.setValue(values["space_group"])
            self._group_changed()
            for key, value in values["lattice"].items():
                self.lattice_controls[key].setValue(value)
            for site in values["sites"]:
                self.add_site(site)
            self.config_path = path
            self._remember(path)
        finally:
            self._updating = False
        self._edited()

    def save_config(self, *_args: Any) -> None:
        try:
            config = self.config()
        except (RuntimeError, TypeError, ValueError) as exc:
            self.failed.emit(str(exc))
            return
        selected, _ = QFileDialog.getSaveFileName(
            self, "Save wyckoff.json", str(self._directory(self.input_path.parent) / "wyckoff.json"),
            "JSON files (*.json)",
        )
        if selected:
            path = Path(selected).expanduser().resolve()
            try:
                path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
            except OSError as exc:
                self.failed.emit(str(exc))
            else:
                self._remember(path)
