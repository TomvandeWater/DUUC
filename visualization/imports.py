import sys
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget,
                             QGridLayout, QVBoxLayout, QHBoxLayout, QLabel, QTextEdit,
                             QLineEdit, QGroupBox, QFormLayout, QComboBox, QMessageBox, QMenu, QSystemTrayIcon)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QIcon, QAction, QColor, QTextCharFormat
from pyvistaqt import QtInteractor
from analysis_modules.aerodynamic import *
from visualize_aircraft import *
from plot_manager import PlotManager
from double_slide import DoubleSlider
import config
import data.atr_reference as ref
import ctypes
from calculation_manager import calculation_manager
import copy
from console_output import *
from analysis_modules.xfoil_run import create_new_polars
