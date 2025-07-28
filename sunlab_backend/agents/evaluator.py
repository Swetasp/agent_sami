from dataclasses import dataclass

@dataclass
class EvaluationResult:
    passed: bool
    explanation: str
    score: float = 1.0  # optional: can default to 1.0 if not using scoring

class Evaluator:
    def __init__(self, llm):
        self.llm = llm

    def evaluate_step(self, step_description, code, stdout, stderr) -> EvaluationResult:
        """
        Check if the step executed successfully and output makes sense.
        This is a basic evaluator; you can later enhance it to use the LLM.
        """
        if "Traceback" in stderr or "Error" in stderr:
            return EvaluationResult(
                passed=False,
                explanation="Error during execution.",
                score=0.0,
            )
        if not stdout.strip():
            return EvaluationResult(
                passed=False,
                explanation="No output produced.",
                score=0.5,
            )
        return EvaluationResult(
            passed=True,
            explanation="Step executed successfully.",
            score=1.0,
        )

