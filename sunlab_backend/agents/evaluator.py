class Evaluator:
    def __init__(self, llm):
        self.llm = llm

    def evaluate(self, code, result):
        """
        Check if output is valid or if errors occurred.
        """
        if "Traceback" in result:
            return False, "Error in execution"
        return True, "Step successful"

