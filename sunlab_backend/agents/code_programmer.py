# agents/code_programmer.py
from dataclasses import dataclass
from typing import Optional
from langchain_core.prompts import PromptTemplate
from utils.json_utils import extract_and_parse_json

@dataclass
class GeneratedCode:
    code: str
    commentary: Optional[str] = None

CODE_GEN_PROMPT = """You are a coding assistant. 
Write **ONLY** valid JSON. Do not add any prose before or after the JSON.

You will receive:
- A step to execute from a biological data analysis pipeline.
- A CSV/H5AD path you can read with scanpy/pandas if needed.

Return a JSON object with this exact schema:

{{
  "code": "A complete python script that can be run as-is",
  "commentary": "Short explanation of what the code does (optional)"
}}

Rules:
- Respond with valid JSON only (no markdown fences).
- Escape all double quotes properly.
- Do not include backticks.

Step to implement: "{step_description}"
Data path: "{data_path}"
"""

class CodeProgrammer:
    def __init__(self, llm):
        self.llm = llm

    def generate(self, step_description: str, data_path: str) -> GeneratedCode:
        prompt = PromptTemplate.from_template(CODE_GEN_PROMPT).format(
            step_description=step_description,
            data_path=data_path,
        )
        raw = self.llm.invoke(prompt)
        # Log raw LLM content for debugging
        print("\n[CodeProgrammer raw LLM output]\n", raw, "\n")

        parsed = extract_and_parse_json(raw)
        print("[CodeProgrammer] Parsed JSON:", parsed)
        if not parsed or "code" not in parsed:
            # Try once to "repair" – or fall back
            repaired = self._best_effort_repair(raw)
            if not repaired:
                raise ValueError("LLM did not return valid JSON for code generation.")
            parsed = repaired

        return GeneratedCode(code=parsed["code"], commentary=parsed.get("commentary"))

    def _best_effort_repair(self, raw: str) -> Optional[dict]:
        """
        Super-light repair: if user forgot JSON, try to wrap the whole thing.
        You can make this smarter (regex, heuristics) if you want.
        """
        raw = raw.strip()
        if raw.startswith("code:"):
            # Example dumb heuristic
            return {"code": raw.split("code:", 1)[1].strip()}
        return None
