from __future__ import annotations

import json
import logging
import re
from typing import Any, Optional

logger = logging.getLogger(__name__)

FENCED_JSON_RE = re.compile(
    r"```(?:json)?\s*(\{.*?\}|\[.*?\])\s*```",
    re.DOTALL | re.IGNORECASE,
)

TOP_LEVEL_JSON_RE = re.compile(
    r"(\{.*\}|\[.*\])",
    re.DOTALL,
)


def extract_and_parse_json(response: str) -> Optional[Any]:
    """
    Try to extract and parse a JSON object/array from an LLM response.

    Strategy (in order):
      1) Look for ```json ... ``` or ``` ... ``` fenced blocks.
      2) Fall back to grabbing the first {...} or [...] top-level block.
      3) As a last resort, slice from the first '{' to the last '}'.

    Returns:
        Parsed JSON (dict/list) on success, or None on failure.
    """
    if not isinstance(response, str):
        logger.warning("LLM response is not a string. Got: %s", type(response))
        return None

    # 1) Try fenced code blocks first
    match = FENCED_JSON_RE.search(response)
    if match:
        json_str = match.group(1).strip()
        parsed = _safe_load(json_str)
        if parsed is not None:
            return parsed

    # 2) Try any top-level JSON-looking block
    match = TOP_LEVEL_JSON_RE.search(response)
    if match:
        json_str = match.group(1).strip()
        parsed = _safe_load(json_str)
        if parsed is not None:
            return parsed

    # 3) Last-resort: slice from first '{' to last '}'
    try:
        start = response.index("{")
        end = response.rindex("}")
        json_str = response[start : end + 1]
        parsed = _safe_load(json_str)
        if parsed is not None:
            return parsed
    except ValueError:
        # '{' or '}' not found
        pass

    logger.error("Failed to extract valid JSON from LLM response.")
    return None


def _safe_load(text: str) -> Optional[Any]:
    """json.loads with logging, returning None on failure."""
    try:
        return json.loads(text)
    except json.JSONDecodeError as e:
        logger.debug("JSON decode failed: %s\nText: %s", e, text[:500])
        return None
