from __future__ import annotations

from typing import Iterable
import xml.etree.ElementTree as ET


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def text(value) -> str:
    if isinstance(value, bool):
        return _bool_text(value)
    return str(value)


def add_text(parent: ET.Element, tag: str, value) -> ET.Element:
    elem = ET.SubElement(parent, tag)
    elem.text = text(value)
    return elem


def add_list_elems(parent: ET.Element, tag: str, values: Iterable) -> ET.Element:
    elem = ET.SubElement(parent, tag)
    for value in values:
        add_text(elem, "elem", value)
    return elem


def add_space_list(parent: ET.Element, tag: str, values: Iterable) -> ET.Element:
    elem = ET.SubElement(parent, tag)
    elem.text = " ".join(text(value) for value in values)
    return elem


def parse_fragment(xml_fragment: str) -> ET.Element:
    return ET.fromstring(xml_fragment)


__all__ = [
    "add_text",
    "add_list_elems",
    "add_space_list",
    "parse_fragment",
    "text",
]
