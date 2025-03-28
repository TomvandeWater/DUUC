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

        if message.startswith("Error:"):
            format.setForeground(QColor("red"))
        else:
            format.setForeground(QColor("black"))

        cursor.insertText(message, format)
        self.text_widget.setTextCursor(cursor)
        self.text_widget.ensureCursorVisible()
        QApplication.processEvents()

    def flush(self):
        pass