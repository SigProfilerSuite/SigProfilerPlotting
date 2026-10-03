"""Interactive SBS-4608 (pentanucleotide + strand context) plotting engine.

Ported from a JS toolkit (AnnotationBuilder.js, BarTraceBuilder.js,
BiasTraceBuilder.js, HeatmapTraceBuilder.js, LayoutBuilder.js,
MutationData.js, SBS4608Data.js, SBS4608PlotBuilder.js,
SBS4608PlotConfig.js, TickLabelFormatter.js, and TraceFormatter.js. The engine
produces plain Plotly-compatible trace/layout dictionaries; the interactive
renderer serializes them to HTML, while the static renderer validates and
exports them through Plotly and Kaleido.
"""
