from PyQt6.QtWidgets import QSlider
from PyQt6.QtCore import pyqtSignal


class DoubleSlider(QSlider):
    doubleValueChanged = pyqtSignal(float)

    def __init__(self, decimals=3, *args, **kwargs):
        super(DoubleSlider, self).__init__(*args, **kwargs)
        self._multiplier = 10 ** decimals
        self.valueChanged.connect(self.emitDoubleValueChanged)

    def emitDoubleValueChanged(self):
        value = float(super(DoubleSlider, self).value()) / self._multiplier
        self.doubleValueChanged.emit(value)

    def value(self):
        return float(super(Slider, self).value()) / self._multiplier

    def setMinimum(self, value):
        return super(DoubleSlider, self).setMinimum(int(value * self._multiplier))

    def setMaximum(self, value):
        return super(DoubleSlider, self).setMaximum(int(value * self._multiplier))

    def setSingleStep(self, value):
        return super(DoubleSlider, self).setSingleStep(int(value * self._multiplier))

    def singleStep(self):
        return float(super(DoubleSlider, self).singleStep()) / self._multiplier

    def setValue(self, value):
        super(DoubleSlider, self).setValue(int(value * self._multiplier))