"""Shared Qt widgets for responsive parameter cards."""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QDoubleSpinBox, QLabel, QSizePolicy


class OneDecimalDisplaySpinBox(QDoubleSpinBox):
    """Display one decimal without rounding the stored calculation value.

    Qt's decimals property also rounds value(), so keep its full double
    precision and format only the text. Untouched display text must never be
    parsed back into the value on focus loss, stepping or config collection.
    """

    def __init__(self) -> None:
        self._text_edited = False
        super().__init__()
        self.setDecimals(323)
        self.setKeyboardTracking(False)
        self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.lineEdit().textEdited.connect(self._mark_edited)
        self.valueChanged.connect(self._finish_editing)
        self.editingFinished.connect(self._finish_editing)

    def textFromValue(self, value: float) -> str:
        # Use the same rounding as XYZ config_type, with the locale's separator.
        rounded = float(f"{value:.1f}")
        return self.locale().toString(rounded, "f", 1)

    def valueFromText(self, text: str) -> float:
        if not self._text_edited:
            return self.value()
        return super().valueFromText(text)

    def setValue(self, value: float) -> None:
        self._text_edited = False
        super().setValue(value)

    def interpretText(self) -> None:
        if self._text_edited:
            super().interpretText()
        self._finish_editing()

    def _mark_edited(self, _text: str) -> None:
        self._text_edited = True

    def _finish_editing(self) -> None:
        self._text_edited = False
        self.lineEdit().setText(
            self.prefix() + self.textFromValue(self.value()) + self.suffix()
        )

    def wheelEvent(self, event: QWheelEvent) -> None:
        event.accept()


class WrappingLabel(QLabel):
    """Let Qt allocate wrapped rows before positioning subsequent controls.

    QLabel already implements heightForWidth(), including font metrics,
    margins and selectable text. Keep that calculation intact: changing
    minimumHeight from resizeEvent can enlarge the label after its form has
    already allocated the next row, causing the two widgets to overlap.
    """

    def minimumSizeHint(self) -> QSize:
        hint = super().minimumSizeHint()
        if self.wordWrap():
            hint.setWidth(0)
        return hint

    def setWordWrap(self, enabled: bool) -> None:
        super().setWordWrap(enabled)
        policy = self.sizePolicy()
        if enabled:
            # Text adapts to the assigned width instead of widening the card.
            policy.setHorizontalPolicy(QSizePolicy.Policy.Ignored)
            policy.setVerticalPolicy(QSizePolicy.Policy.Preferred)
            policy.setHeightForWidth(True)
        else:
            policy.setHorizontalPolicy(QSizePolicy.Policy.Preferred)
            policy.setVerticalPolicy(QSizePolicy.Policy.Preferred)
            policy.setHeightForWidth(False)
        self.setSizePolicy(policy)
        self.updateGeometry()
