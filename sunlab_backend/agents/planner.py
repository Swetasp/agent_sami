class Planner:
    def __init__(self, llm):
        self.llm = llm

    def plan(self, task: str, data_signature: str):
        """
        Convert a high-level task into a list of steps.
        """
        prompt = f"Task: {task}\nData: {data_signature}\nBreak this into steps."
        response = self.llm.invoke(prompt)
        # For now, return a dummy plan
        return [
            {"id": 1, "description": "Normalize data"},
            {"id": 2, "description": "Perform clustering"},
        ]
