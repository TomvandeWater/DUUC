import sys
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget, QPushButton, QFileDialog,
                             QGridLayout, QVBoxLayout, QHBoxLayout, QLabel, QTextEdit, QInputDialog,
                             QLineEdit, QGroupBox, QFormLayout, QComboBox, QMessageBox, QMenu, QSystemTrayIcon)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QIcon, QAction, QColor, QTextCharFormat, QPixmap, QImage, QFont, QDoubleValidator
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
import json
import yaml
from termcolor import colored
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


def convert_for_json(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert numpy arrays to lists
    elif isinstance(obj, dict):
        return {k: convert_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_for_json(i) for i in obj]
    else:
        return obj