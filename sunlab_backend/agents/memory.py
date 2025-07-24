class GlobalMemory:
    """
// Global memory stores the final code for each historical step.
    """
    def __init__(self):
        self.codes = []

    def add_code(self, code):
        self.codes.append(code)

    def get_all_code(self):
        return "\n".join(self.codes)
