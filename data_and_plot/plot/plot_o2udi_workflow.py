#!/usr/bin/env python3
"""Generate the O2UDI overview as a PNG with the Python Graphviz library.

Requirements: Python packages ``graphviz`` and ``cairosvg``, and a Graphviz
installation that includes ``neato``. Install the Python packages with
``pip install graphviz cairosvg``. CairoSVG also requires the Cairo library.

Usage:
    python plot_o2udi_workflow.py
    python plot_o2udi_workflow.py --output-dir figures --dpi 400 --width-mm 180

The default output directory is the script directory. The only output file is
fig_o2udi_workflow.png. Graphviz and SVG data remain in memory. The diagram
uses short block labels, two pastel module groups, and one inference loop.
Positions are in inches and can be adjusted in ``make_graph``.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess

from graphviz import Digraph

STEM = "fig_o2udi_workflow"
# Pango treats a bare trailing "Roman" as a style. The comma terminates the
# family name, so Windows Graphviz does not interpret it as "Times New".
FONT = "Times New Roman,"
INK = "#252525"
GROUP_LINE = "#425C80"
BLUE = "#BBCDE9"
PEACH = "#FCE7DA"
ORANGE = "#F6CCAE"
PALE_YELLOW = "#FFF3CC"
YELLOW = "#FFE79C"


def make_graph(dpi: int) -> Digraph:
    """Draw a compact overview; detailed equations remain in the manuscript.

    The unobservable-state block is initialized before inference and receives
    the updated state on each cycle. The historical window is refreshed at
    the new epoch. Scaling transfers the measured mean orbital-change rates;
    mean perigee precession supplies the complementary angular update.
    """
    g = Digraph("O2UDI_workflow", engine="neato")
    g.attr(
        "graph",
        bgcolor="white",
        fontname=FONT,
        pad="0.14",
        margin="0",
        overlap="true",
        splines="polyline",
        outputorder="edgesfirst",
        dpi=str(dpi),
    )
    g.attr(
        "node",
        shape="box",
        style="rounded,filled",
        fontname=FONT,
        fontsize="18",
        fontcolor=INK,
        color=INK,
        penwidth="0.9",
        margin="0.15,0.10",
        width="2.65",
        height="0.60",
        pin="true",
    )
    g.attr(
        "edge",
        color=INK,
        penwidth="1.0",
        arrowsize="0.68",
        fontname=FONT,
        fontsize="17",
    )

    def block(graph: Digraph,
              name: str,
              label: str,
              x: float,
              y: float,
              fill: str,
              *,
              width: float = 2.65,
              height: float = 0.60) -> None:
        graph.node(name,
                   label,
                   pos=f"{x},{y}!",
                   fillcolor=fill,
                   width=str(width),
                   height=str(height))

    def point(graph: Digraph, name: str, x: float, y: float) -> None:
        graph.node(name,
                   label="",
                   shape="point",
                   style="invis",
                   width="0.001",
                   height="0.001",
                   pos=f"{x},{y}!")

    def words(graph: Digraph,
              name: str,
              label: str,
              x: float,
              y: float,
              *,
              bold: bool = False,
              size: int = 17) -> None:
        if bold:
            label = "<<B>" + label.replace("\n", "<BR/>") + "</B>>"
        graph.node(name,
                   label=label,
                   shape="plain",
                   style="",
                   width="0",
                   height="0",
                   margin="0",
                   fontsize=str(size),
                   pos=f"{x},{y}!")

    # Cluster corner anchors reserve padding explicitly: neato does not reserve
    # dot-style cluster-label space around pinned nodes. Titles are separate
    # short text nodes, positioned away from all incoming and outgoing arrows.
    with g.subgraph(name="cluster_matching") as c:
        c.attr(label="",
               style="rounded,filled,dashed",
               color=GROUP_LINE,
               fillcolor=PEACH,
               penwidth="0.9")
        point(c, "match_corner_sw", 0.65, 4.20)
        point(c, "match_corner_ne", 3.85, 6.95)
        words(c,
              "match_title",
              "Dynamic fragment\nmatching",
              2.25,
              6.54,
              bold=True,
              size=21)
        block(c, "radial", "Perigee–apogee matching", 2.25, 5.80, ORANGE)
        block(c, "angular", "Perigee-angle matching", 2.25, 4.85, ORANGE)
        c.edge("radial:s", "angular:n")

    with g.subgraph(name="cluster_inference") as c:
        c.attr(label="",
               style="rounded,filled,dashed",
               color=GROUP_LINE,
               fillcolor=PALE_YELLOW,
               penwidth="0.9")
        point(c, "infer_corner_sw", 0.65, 1.12)
        point(c, "infer_corner_ne", 7.85, 3.87)
        words(c,
              "infer_title",
              "Cross-size orbital inference",
              3.05,
              1.48,
              bold=True,
              size=21)
        block(c,
              "response",
              "Observed-response\nestimation",
              2.25,
              3.24,
              YELLOW,
              height=0.70)
        block(c,
              "scaling",
              "Area-to-mass\nratio scaling",
              2.25,
              2.20,
              YELLOW,
              height=0.70)
        block(c,
              "precession",
              "Mean perigee\nprecession",
              6.15,
              3.24,
              YELLOW,
              width=2.95,
              height=0.70)
        block(c,
              "update",
              "Orbital-state update",
              6.15,
              2.20,
              YELLOW,
              width=2.95)
        c.edge("response:s", "scaling:n")
        c.edge("scaling:e", "update:w")
        c.edge("precession:s", "update:n")

    block(g,
          "observed",
          "Observable-fragment\norbital records",
          4.20,
          7.70,
          BLUE,
          width=3.10,
          height=0.75)
    block(g,
          "target",
          "Unobservable-fragment\nstate",
          6.15,
          5.80,
          BLUE,
          width=2.95,
          height=0.75)
    block(g,
          "output",
          "Inferred orbital histories",
          6.15,
          0.50,
          BLUE,
          width=2.95,
          height=0.62)

    # Observable candidates are restricted to the current historical window.
    # Route this input around, rather than through, the matching-group title.
    point(g, "records_left_top", 0.20, 7.70)
    point(g, "records_left_bottom", 0.20, 5.80)
    g.edge("observed:w", "records_left_top", arrowhead="none")
    g.edge("records_left_top", "records_left_bottom", arrowhead="none")
    g.edge("records_left_bottom", "radial:w")
    words(g, "window_label", "Historical window", 1.40, 7.90)

    g.edge("target:w", "radial:e")
    g.edge("target:s", "precession:n")
    g.edge("angular:s", "response:n")
    g.edge("update:s", "output:n")

    # The updated target state and epoch drive the next matching cycle.
    point(g, "loop_right_bottom", 8.30, 2.20)
    point(g, "loop_right_top", 8.30, 5.80)
    g.edge("update:e", "loop_right_bottom", arrowhead="none")
    g.edge("loop_right_bottom", "loop_right_top", arrowhead="none")
    g.edge("loop_right_top", "target:e")
    words(g, "loop_label", "Next epoch", 8.95, 3.95)
    return g


def graphviz_svg(graph: Digraph) -> str:
    """Render with desktop Graphviz or the available bundled Graphviz/WASM."""
    if shutil.which("neato"):
        return graph.pipe(format="svg", encoding="utf-8")
    # This optional fallback makes the script runnable in the authoring runtime.
    # A normal desktop installation uses neato above and does not need Node.js.
    module_root = os.environ.get("CODEX_PRIMARY_RUNTIME_NODE_MODULES")
    viz_path = Path(module_root or ".") / "@viz-js/viz/dist/viz.js"
    if not module_root or not viz_path.is_file() or not shutil.which("node"):
        raise RuntimeError(
            "Install Graphviz and add its bin directory (containing neato) to PATH."
        )
    bridge = r'''
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
const { instance } = await import(pathToFileURL(process.argv[1]).href);
const viz = await instance();
const result = viz.render(readFileSync(0, "utf8"), { engine: "neato", format: "svg" });
if (result.status !== "success") {
  console.error(JSON.stringify(result.errors));
  process.exit(1);
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
    return result.stdout


def normalize_svg(svg: str, width_mm: float) -> str:
    """Set the in-memory canvas width before PNG rasterization."""
    box = [
        float(v) for v in re.search(r'viewBox="([\d.\s]+)"', svg)[1].split()
    ]
    height_mm = width_mm * box[3] / box[2]
    svg = re.sub(
        r'<svg width="[^"]+" height="[^"]+"',
        f'<svg width="{width_mm:g}mm" height="{height_mm:.4f}mm"',
        svg,
    )
    # The comma above is for Pango's font parser. Use a normal CSS fallback
    # list when CairoSVG rasterizes the result (no SVG file is written).
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
    args.output_dir.mkdir(parents=True, exist_ok=True)
    graph = make_graph(args.dpi)
    import cairosvg

    svg = normalize_svg(graphviz_svg(graph), args.width_mm)
    png_path = args.output_dir / f"{STEM}.png"
    cairosvg.svg2png(bytestring=svg.encode(),
                     dpi=args.dpi,
                     write_to=str(png_path))
    print(f"Saved PNG: {png_path.resolve()}")


if __name__ == "__main__":
    main()
