# src/planner.py

from typing import List, Dict, Any
from langchain_core.prompts import PromptTemplate
from utils.json_utils import extract_and_parse_json


class Planner:
    """
    The Planner class is responsible for breaking down a high-level user task
    into smaller, executable steps for data analysis.
    """

    def __init__(self, llm, data_representation: str):
        """
        Initialize the Planner.

        Parameters
        ----------
        llm : object
            The large language model instance (e.g., Groq/OpenAI/Ollama via LangChain).
        data_representation : str
            A string description of the dataset (shape, observations, features, source, etc.).
        """
        self.llm = llm
        self.data_representation = data_representation

        # Prompt template used for task planning
        self.prompt_template = """You are an assistant specialized in bioinformatics.
Based on the information below, generate a **detailed and executable** task plan for the user's multi-omics or single-cell data analysis.
The output **must** be valid JSON (no extra text), so it can be directly parsed for subtask extraction and execution.

User Task:
{user_task}

Dataset Description:
{data_representation}

Strictly return a JSON object in the following format (no extra fields):
{{
  "steps": [
    {{
      "id": 1,
      "description": "Description of step 1 (clear and concise, suitable for a code generator)."
    }},
    {{
      "id": 2,
      "description": "Description of step 2"
    }}
  ]
}}
"""

    def plan(self, user_task: str) -> List[Dict[str, Any]]:
        """
        Break down the high-level task into a list of executable steps.

        Returns
        -------
        List[Dict[str, Any]]
            A list of steps, e.g., [{"id": 1, "description": "..."}, ...].
        """
        prompt = PromptTemplate(
            input_variables=["user_task", "data_representation"],
            template=self.prompt_template,
        ).format(
            user_task=user_task,
            data_representation=self.data_representation,
        )

        # Use the LLM to get the response
        response = self.llm.invoke(prompt)
        print("Planner response - " + response)
        # Attempt to extract and parse the JSON
        parsed_json = extract_and_parse_json(response)

        if parsed_json and isinstance(parsed_json, dict):
            steps = parsed_json.get("steps", [])
            if isinstance(steps, list):
                return steps

        # Fallback: if parsing fails, return a default plan
        print("Failed to parse the planner output as JSON, returning default steps.")
        return [
            {"id": 1, "description": "Normalize data"},
            {"id": 2, "description": "Perform clustering"},
        ]

    def generate_final_result(self):
        """
        Generate the final result summary.
        Since all code and results are stored in a Jupyter Notebook, 
        this can simply return the notebook path or a summary of the analysis.
        """
        pass
