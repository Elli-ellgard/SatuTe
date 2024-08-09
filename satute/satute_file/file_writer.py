class FileWriter:
    def __init__(self, filename):
        self.filename = filename
        self.file = None

    def open_file(self, mode="w"):
        self.file = open(self.filename, mode)

    def write_to_file(self, content):
        if self.file:
            self.file.write(content)
        else:
            raise ValueError("File is not open. Call open_file() before writing.")

    def write_header(self, content):
        if self.file:
            self.file.write(f"\n#### {content} ####\n")
        else:
            raise ValueError("File is not open. Call open_file() before writing.")

    def close_file(self):
        if self.file:
            self.file.close()
            self.file = None

    def __enter__(self):
        self.open_file()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close_file()
