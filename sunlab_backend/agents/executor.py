# agents/executor.py
import tempfile
import subprocess
import textwrap
from typing import Optional
from agents.code_programmer import CodeProgrammer, GeneratedCode
from utils.json_utils import extract_and_parse_json
from code_sandbox import CodeSandbox  # your working sandbox
from dataclasses import dataclass
@dataclass
class ExecutionResult:
    success: bool
    stdout: str
    stderr: str
    notebook_path: Optional[str]
    code: Optional[str] = None
    # explanation: str
class Executor:
    def __init__(self, llm, data_path: str, memory, sandbox: CodeSandbox):
        self.llm = llm
        self.data_path = data_path
        self.memory = memory
        self.programmer = CodeProgrammer(llm)
        self.sandbox = sandbox
    def _maybe_pip_install(self, requirements):
        if not requirements:
            return
        cmd = ["python", "-m", "pip", "install", *requirements]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            # Non-fatal: we’ll still try to run the code
            print("[Executor] pip install failed:", e.stderr)
    def _run_python_locally(self, code: str, timeout: int = 600):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
            f.write(code)
            path = f.name
        proc = subprocess.Popen(
            ["python", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
            success = proc.returncode == 0
        except subprocess.TimeoutExpired:
            proc.kill()
            return False, "", "TimeoutExpired", path
        return success, stdout, stderr, path
    def execute_step(self, step_description: str, user_requirements, data_description, global_memory) -> ExecutionResult:
        # 1) Ask LLM for code
        gen: GeneratedCode = self.programmer.generate(step_description, self.data_path, user_requirements, data_description, global_memory)
        # 2) Install optional requirements (best effort)
        # self._maybe_pip_install(gen.requirements)
        # 3) (Optional) also drop the code cell in a notebook
        self.sandbox.add_markdown_cell(f"### Step: {step_description}\n")
        self.sandbox.add_code_cell(gen)
        # 4) Execute locally (fast feedback)
        success, stdout, stderr, _ = self._run_python_locally(gen)
        nb_path = None
        # 5) Keep a record in the notebook (so user can open / replay)
        #    You can decide to only save if success, or always save.
        # 
        # 6) Persist result into memory
        # self.memory.add({"step": step_description, "success": success, "stdout": stdout, "stderr": stderr})
        return ExecutionResult(
            success=success,
            stdout=stdout,
            stderr=stderr,
            notebook_path=nb_path,
            code=gen
            # explanation=gen.explanation
        )
    def finalize_notebook(self) -> Optional[str]:
            try:
                result = self.sandbox.execute()
                return result.notebook_path
            except Exception as e:
                print("[Executor] Notebook execution failed:", e)
                return None


















