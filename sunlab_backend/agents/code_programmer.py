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
Please write Python code to complete the current step based on the following information.

You will receive:
- A step to execute from a biological data analysis pipeline.
- A CSV/H5AD path you can read with scanpy/pandas if needed.

Return a JSON object with this exact schema:

Step to implement: "{step_description}"
Data path: "{data_path}"

Rules
- The code must use the data file path: {data_path}
- Always use pandas to load CSV files, and use anndata.AnnData(df) to convert to AnnData if needed.
- The code should include necessary imports and data loading.
- Do not add any comments, just provide the code.
- The output format should be a code block enclosed in ```python and ```
- ALWAYS treat file paths as raw strings (prefix them with r"") to avoid invalid escape sequence errors on Windows.
- Use only functions, classes, or methods that are available and recommended in the latest stable version of any library used (e.g., scanpy, pandas, numpy, etc.).
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
        # print("RAW: ")
        # print(raw)
        code = self.extract_code(raw.content)
        # Log raw LLM content for debugging
        print("\n[CodeProgrammer raw LLM output]\n", raw.content, "\n")
        print("\n[CodeProgrammer code LLM output]\n", code, "\n")

        # parsed = extract_and_parse_json(raw)
        # print("[CodeProgrammer] Parsed JSON:", parsed)
        # if not parsed or "code" not in parsed:
        #     # Try once to "repair" – or fall back
        #     repaired = self._best_effort_repair(raw)
        #     if not repaired:
        #         raise ValueError("LLM did not return valid JSON for code generation.")
        #     parsed = repaired

        # return GeneratedCode(code=parsed["code"], commentary=parsed.get("commentary"))
        return code

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

    def extract_code(self, response):
        start = response.find('```python')
        end = response.find('```', start + 9)
        if start != -1 and end != -1:
            return response[start+9:end].strip()
        else:
            return response.strip()