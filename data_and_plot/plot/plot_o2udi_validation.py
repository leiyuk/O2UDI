#!/usr/bin/env python3
"""Draw the O2UDI experimental-validation overview as one PNG.

Requirements:
    pip install graphviz cairosvg
    Install Graphviz and put its bin directory (containing dot) on PATH.
    CairoSVG also requires Cairo, as in plot_o2udi_workflow.py.

Usage:
    python plot_o2udi_validation.py
    python plot_o2udi_validation.py --output-dir figures --dpi 400 --width-mm 180

Only fig_o2udi_validation.png is written. DOT and SVG remain in memory.
The dot engine lays out three side-by-side validation paths with a shared
evidence box below them; no neato triangulation or pinned coordinates are needed.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess

from graphviz import Digraph

STEM = "fig_o2udi_validation"
# A terminating comma prevents Pango from reading "Roman" as a font style.
FONT = "Times New Roman,"
INK = "#252525"
GROUP_LINE = "#425C80"
EVIDENCE_LABEL = "Complementary evidence for O2UDI"

# Edit the labels here. Only the two comparisons use three lines, with vs.
# on its own line; the regression description stays on one line.
# Titles, inputs, and validation scopes stay on one line; no LaTeX is required.
# These are parallel sources of evidence, not sequential algorithm stages.
EXPERIMENTS = (
    {
        "key": "numerical",
        "title": "Closed numerical experiment",
        "background": "#FCE7DA",
        "fill": "#F6CCAE",
        "input": "Simulated large-fragment orbits as O2UDI input",
        "comparison":
        "Inferred vs. numerical reference\npopulations of small fragments",
        "scope": "Distribution accuracy and orbital structure",
    },
    {
        "key": "cosmos1408",
        "title": "Cosmos 1408 TLE analysis",
        "background": "#E9F0FA",
        "fill": "#BBCDE9",
        "input": "Mean semimajor-axis decay rate and B*",
        "comparison": "Linear regression for each fragment",
        "scope": "Fragment-level response proportionality",
    },
    {
        "key": "cosmos2251",
        "title": "Cosmos 2251 TLE experiment",
        "background": "#FFF3CC",
        "fill": "#FFE79C",
        "input": "Lower-mean-B* group as O2UDI input",
        "comparison":
        "Inferred vs. true\ndistributions of the higher-mean-B* group",
        "scope": "Long-term cross-size population inference",
    },
)


def bold_label(label: str) -> str:
    """Use Graphviz HTML only for bold text and explicit line breaks."""
    from html import escape

    return "<<B>" + escape(label).replace("\n", "<BR/>") + "</B>>"


def make_graph(dpi: int) -> Digraph:
    """Arrange three panels horizontally above a shared evidence box.

    The numerical reference is obtained by direct propagation, not by O2UDI.
    Cosmos 1408 supports the fragment-level scaling relation through TLE
    regression; it is not another population-inference run.
    Both Cosmos 2251 groups are cataloged. The higher-mean-B* group's orbital
    history is withheld from inference and supplies the true distribution.
    """
    graph = Digraph("O2UDI_validation", engine="dot")
    graph.attr(
        "graph",
        bgcolor="white",
        rankdir="TB",
        newrank="true",
        splines="ortho",
        nodesep="0.90",
        ranksep="0.50",
        pad="0.14",
        margin="0",
        outputorder="edgesfirst",
        fontname=FONT,
        dpi=str(dpi),
    )
    graph.attr(
        "node",
        shape="box",
        style="rounded,filled",
        fontname=FONT,
        fontsize="24",
        fontcolor=INK,
        color=INK,
        penwidth="0.9",
        width="5.60",
        height="0.60",
        margin="0.16,0.12",
    )
    graph.attr("edge", color=INK, penwidth="1.0", arrowsize="0.70")

    for experiment in EXPERIMENTS:
        key = experiment["key"]
        with graph.subgraph(name=f"cluster_{key}") as panel:
            panel.attr(
                label=bold_label(experiment["title"]),
                fontname=FONT,
                fontsize="28",
                fontcolor=INK,
                labelloc="t",
                labeljust="c",
                style="rounded,filled,dashed",
                color=GROUP_LINE,
                fillcolor=experiment["background"],
                penwidth="0.9",
                margin="20",
            )
            panel.node(f"{key}_input",
                       experiment["input"],
                       fillcolor=experiment["fill"],
                       group=key)
            panel.node(f"{key}_comparison",
                       experiment["comparison"],
                       fillcolor=experiment["fill"],
                       group=key,
                       height="0.90")
            panel.node(f"{key}_scope",
                       bold_label(experiment["scope"]),
                       fillcolor=experiment["fill"],
                       group=key,
                       height="0.62")
            panel.edge(f"{key}_input:s", f"{key}_comparison:n")
            panel.edge(f"{key}_comparison:s", f"{key}_scope:n")

    # Align matching rows without changing any panel labels or node styles.
    for stage in ("input", "comparison", "scope"):
        with graph.subgraph(name=f"rank_{stage}") as row:
            row.attr(rank="same")
            for experiment in EXPERIMENTS:
                row.node(f"{experiment['key']}_{stage}")

    # Invisible edges preserve the left-to-right order of the three panels.
    for left, right in zip(EXPERIMENTS, EXPERIMENTS[1:]):
        graph.edge(f"{left['key']}_input",
                   f"{right['key']}_input",
                   style="invis",
                   weight="100",
                   constraint="false",
                   arrowhead="none")

    graph.node("evidence",
               bold_label(EVIDENCE_LABEL),
               fillcolor="#F2F4F7",
               fontsize="28",
               width="22.0",
               height="0.72")
    for experiment in EXPERIMENTS:
        graph.edge(f"{experiment['key']}_scope:s",
                   "evidence",
                   minlen="2")
    return graph


def graphviz_svg(graph: Digraph) -> str:
    """Render in memory using desktop Graphviz, with an authoring fallback."""
    if shutil.which("dot"):
        return graph.pipe(format="svg", encoding="utf-8")

    # Optional fallback in the authoring environment. A normal desktop uses
    # Graphviz above and does not require Node.js or the WASM package.
    module_root = os.environ.get("CODEX_PRIMARY_RUNTIME_NODE_MODULES")
    viz_path = Path(module_root or ".") / "@viz-js/viz/dist/viz.js"
    if not module_root or not viz_path.is_file() or not shutil.which("node"):
        raise RuntimeError(
            "Install Graphviz and add its bin directory (containing dot) to PATH."
        )
    bridge = r'''
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
const { instance } = await import(pathToFileURL(process.argv[1]).href);
const viz = await instance();
const result = viz.render(readFileSync(0, "utf8"), {engine: "dot", format: "svg"});
if (result.status !== "success") {
  console.error(JSON.stringify(result.errors));
  process.exit(1);
}
for (const message of result.errors ?? []) {
  console.error(message.message);
}
process.stdout.write(result.output);
'''
    result = subprocess.run(
        ["node", "--input-type=module", "-e", bridge,
         str(viz_path)],
        input=graph.source,
        text=True,
        capture_output=True,
        check=True,
    )
    if result.stderr:
        import sys

        print(result.stderr.rstrip(), file=sys.stderr)
    return result.stdout


def normalize_svg(svg: str, width_mm: float) -> str:
    """Set the final physical width and a regular CSS font fallback list."""
    match = re.search(r'viewBox="([\d.\s+-]+)"', svg)
    if match is None:
        raise ValueError("Graphviz returned an SVG without a viewBox.")
    box = [float(value) for value in match[1].split()]
    height_mm = width_mm * box[3] / box[2]
    svg = re.sub(
        r'<svg width="[^"]+" height="[^"]+"',
        f'<svg width="{width_mm:g}mm" height="{height_mm:.4f}mm"',
        svg,
        count=1,
    )
    return svg.replace(
        'font-family="Times New Roman,"',
        'font-family="Times New Roman, Nimbus Roman, Liberation Serif, serif"',
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",
                        type=Path,
                        default=Path(__file__).resolve().parent)
    parser.add_argument("--dpi", type=int, default=400)
    parser.add_argument("--width-mm", type=float, default=180)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.width_mm <= 0:
        parser.error("--width-mm must be positive")

    import cairosvg

    svg = normalize_svg(graphviz_svg(make_graph(args.dpi)), args.width_mm)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"{STEM}.png"
    cairosvg.svg2png(bytestring=svg.encode("utf-8"),
                     dpi=args.dpi,
                     background_color="#FFFFFF",
                     write_to=str(output))
    print(f"Saved PNG: {output.resolve()}")


if __name__ == "__main__":
    main()
