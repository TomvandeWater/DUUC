from imports import *


class ConsoleOutput:
    def __init__(self, text_widget):
        self.text_widget = text_widget
        self.stdout = sys.stdout
        self.stderr = sys.stderr

    def write(self, message):
        if self.text_widget.textCursor().atEnd():
            self.text_widget.moveCursor(self.text_widget.textCursor().End)

        cursor = self.text_widget.textCursor()
        format = QTextCharFormat()

        if "Matlab" in message:
            format.setForeground(QColor(Qt.blue))  # Red for "Matlab"
        elif "Error" in message:
            format.setForeground(QColor(Qt.red))  # Red for error messages
        elif "Success" in message:
            format.setForeground(QColor(Qt.darkGreen))  # Green for success messages
        elif "XFoil" in message:
            format.setForeground(QColor(Qt.magenta))
        elif "plots" in message:
            format.setForeground(QColor(Qt.darkGray))

        cursor.insertText(message, format)
        self.text_widget.setTextCursor(cursor)
        self.text_widget.ensureCursorVisible()
        QApplication.processEvents()

    def flush(self):
        pass