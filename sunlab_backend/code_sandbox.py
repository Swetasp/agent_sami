from __future__ import annotations

import os
import tempfile
import nbformat
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional, Dict, Any

from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError


@dataclass
class ExecutionResult:
    success: bool
    notebook_path: str
    stdout: str = ""
    stderr: str = ""
    errors: List[str] = field(default_factory=list)
    cell_outputs: List[Dict[str, Any]] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.success


class CodeSandbox:
    """
    A lightweight Jupyter-Notebook-based code sandbox to run dynamically generated
    analysis code safely and persist the full execution trace to disk.
    """

    def __init__(
        self,
        notebook_path: str,
        kernel_name: str = "python3",
        timeout: int = 1200,
        allow_errors: bool = False,
        workdir: Optional[str] = None,
    ):
        """
        Parameters
        ----------
        notebook_path : str
            Where to save the executed notebook ('.ipynb').
        kernel_name : str
            The Jupyter kernel to use.
        timeout : int
            Per-cell timeout (seconds).
        allow_errors : bool
            If True, the notebook continues executing after errors.
        workdir : Optional[str]
            Directory to execute the notebook in. Defaults to a temp dir.
        """
        self.notebook_path = notebook_path
        self.kernel_name = kernel_name
        self.timeout = timeout
        self.allow_errors = allow_errors
        self.workdir = workdir or tempfile.mkdtemp(prefix="sunlab_exec_")

        self.nb = nbformat.v4.new_notebook()
        self.nb.metadata["kernelspec"] = {"name": kernel_name, "language": "python"}
        self.nb["cells"] = []

    # ------------- builders -------------------------------------------------

    def add_code_cell(self, code: str) -> None:
        """Append a code cell to the notebook."""
        self.nb["cells"].append(nbformat.v4.new_code_cell(code))

    def add_markdown_cell(self, text: str) -> None:
        """Append a markdown cell to the notebook."""
        self.nb["cells"].append(nbformat.v4.new_markdown_cell(text))

    def add_parameters_cell(self, params: Dict[str, Any]) -> None:
        """
        Inject a parameters cell (Papermill style, but without requiring papermill).
        The dict is turned into Python assignments at the top of the notebook.
        """
        assignments = []
        for k, v in params.items():
            assignments.append(f"{k} = {repr(v)}")
        code = "# Parameters\n" + "\n".join(assignments)
        self.nb["cells"].insert(0, nbformat.v4.new_code_cell(code))

    # ------------- execution ------------------------------------------------

    def execute(self, parameters: Optional[Dict[str, Any]] = None) -> ExecutionResult:
        """Run the notebook and persist it to disk."""
        # Optionally inject parameters
        if parameters:
            self.add_parameters_cell(parameters)

        unique_notebook_path = self._unique_path(self.notebook_path)

        ep = ExecutePreprocessor(
            timeout=self.timeout,
            kernel_name=self.kernel_name,
            allow_errors=self.allow_errors,
            store_widget_state=True,
        )

        # Where the notebook will run
        resources = {"metadata": {"path": self.workdir}}

        stdout, stderr = [], []
        errors: List[str] = []
        cell_outputs: List[Dict[str, Any]] = []

        try:
            ep.preprocess(self.nb, resources)

            # Collect outputs
            for cell in self.nb.cells:
                if cell.cell_type != "code":
                    continue
                co = {"source": cell.source, "outputs": []}
                for out in cell.get("outputs", []):
                    co["outputs"].append(out)
                    if out.get("name") == "stdout" and "text" in out:
                        stdout.append(out["text"])
                    if out.get("name") == "stderr" and "text" in out:
                        stderr.append(out["text"])
                    if out.get("output_type") == "error":
                        errors.append(
                            f"{out.get('ename', '')}: {out.get('evalue', '')}\n{''.join(out.get('traceback', []))}"
                        )
                cell_outputs.append(co)

        except CellExecutionError as e:
            # Save whatever executed so far
            errors.append(str(e))
        finally:
            os.makedirs(os.path.dirname(unique_notebook_path) or ".", exist_ok=True)
            with open(unique_notebook_path, "w", encoding="utf-8") as f:
                nbformat.write(self.nb, f)

        return ExecutionResult(
            success=len(errors) == 0,
            notebook_path=unique_notebook_path,
            stdout="".join(stdout),
            stderr="".join(stderr),
            errors=errors,
            cell_outputs=cell_outputs,
        )

    # ------------- utils ----------------------------------------------------

    @staticmethod
    def _unique_path(path: str) -> str:
        """Add a timestamp if the path already exists."""
        if not os.path.exists(path):
            return path
        base, ext = os.path.splitext(path)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"{base}_{ts}{ext}"
if __name__ == "__main__":
    sandbox = CodeSandbox("/Users/swetapasari/Downloads/SUN Lab/analysis.ipynb")
    sandbox.add_markdown_cell("# SunLab Auto Analysis")
    sandbox.add_code_cell("import numpy as np\nprint('Hello from sandbox!', np.arange(5))")
    result = sandbox.execute(parameters={"alpha": 0.1, "dataset": "brain_glycomics.csv"})
    print("OK:", result.success, "Saved to:", result.notebook_path)
