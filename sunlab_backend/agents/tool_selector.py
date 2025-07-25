import json

class ToolSelector:
    """
    Asks the LLM which tool and parameters should be used for a given step.
    """
    def __init__(self, llm, tools):
        self.llm = llm
        self.tools = tools  # List of {name, description}

    def select_tool(self, step_description: str):
        """
        Ask the LLM to pick a tool for the given step description.
        """
        tools_list = "\n".join(
            [f"- {tool['name']}: {tool['description']}" for tool in self.tools]
        )

        prompt = f"""
        You are a tool selection agent.
        Task: {step_description}

        Tools available:
        {tools_list}

        Respond with a JSON like:
        {{
          "tool_name": "<tool_name>",
          "parameters": {{}}
        }}
        """

        response = self.llm.invoke(prompt)
        try:
            start = response.index("{")
            end = response.rindex("}") + 1
            return json.loads(response[start:end])
        except Exception as e:
            print(f"[ToolSelector] Failed to parse tool selection: {e}")
            return {"tool_name": "none", "parameters": {}}
