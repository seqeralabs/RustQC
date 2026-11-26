//! Qualimap-compatible HTML report generator.
//!
//! Produces a self-contained `qualimapReport.html` with embedded CSS and image
//! assets, matching Qualimap's Sphinx agogo theme layout exactly.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use super::format::{format_decimal, format_long, format_percentage, percentage_of};
use super::output::{NUM_PRIME_BASES, NUM_TRANSCRIPTS_FOR_BIAS};
use crate::Strandedness;
use log::debug;

// ===================================================================
// Embedded assets (base64-encoded)
// ===================================================================

/// Qualimap logo (small) — base64-encoded PNG.
const LOGO_B64: &str = "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";

/// Header background texture — base64-encoded PNG.
const BGTOP_B64: &str = "iVBORw0KGgoAAAANSUhEUgAAAAMAAABQCAYAAADcMs/mAAAAAXNSR0IArs4c6QAAAAZiS0dEAFUAVwBTkiHDvwAAAAlwSFlzAAALEwAACxMBAJqcGAAAAAd0SU1FB9kCEwctBZSG9LMAAAAZdEVYdENvbW1lbnQAQ3JlYXRlZCB3aXRoIEdJTVBXgQ4XAAABCUlEQVQ4y12Sy45AIQhDD8Rfn2QyH21nQVG8rhR5lLbx+/cjfJJxEoQASeSuGwAL1HcSCfxYFa3Aki8RkNr1kGDhtq5xjmQEz4/PqmYVSD1D+yZIOU0NFEP6dhMi3O3AUbUOB1b97DG0xlQDAOIAjV67k86cystGI4lsdlTYqiYOatMwUEeneZ9WIS42sZv4YEiiU3NWiGa02eECHUP1puX2HPvgynjYQe7Gji8H6ppdhBSz4U15jFSbhdqW28TfOR928vJuteOIZaqqwQaFYGMVHlNMRlVzrem1JYTpSrEnvfZotAp2fL7ubek1a/R1/Hzk8PiXKmw+xCqaNkSM5aS7zxl6ZJzd/gFE3hujlFg/uwAAAABJRU5ErkJggg==";

/// Footer background texture — base64-encoded PNG.
const BGFOOTER_B64: &str = "iVBORw0KGgoAAAANSUhEUgAAAAMAAABQCAYAAADcMs/mAAAAAXNSR0IArs4c6QAAAAZiS0dEANMA1wDPtAtCYgAAAAlwSFlzAAALEwAACxMBAJqcGAAAAAd0SU1FB9kCEwc0ERVcidYAAAAZdEVYdENvbW1lbnQAQ3JlYXRlZCB3aXRoIEdJTVBXgQ4XAAABDUlEQVQ4y22SO5IlMQgEMxXc/2prjDHHYQ0+0nsxbalDoiiS8uf3X9Lf4fk+f+yDWjcKkAQKCQgHBKfGeqetZgJIgFt2BNJ8n5VCVA/rGTptOOY46JsyR0sjkBzsArwjwD5r1445IdK+SzhzKGln7uSMSVMCs2psO67RYeUSzYfBZV0VaY/d4ncLDtHkCyJIjFajqnkEYqDn0MEa9qy4XwuOXGkIt/8TCoWjHyNsKoqb5OWWiNkODpD+HbFqHKPkOmiR3XaNMHHRy9plUHpEPdm4QOagQs5k9CZk+nTZM9ywznzHnoiNU7vPs8Y2CZwG8+S6kxWfAg/4cHf/Rtkv8E+q3uxc6dKI6+xrC/8BSHkzA/8EsJQAAAAASUVORK5CYII=";

// ===================================================================
// Embedded CSS (inlined from agogo.css + report.css)
// ===================================================================

/// Combined CSS from Qualimap's agogo theme + report.css, faithfully reproduced.
/// Background image references replaced with data URIs at render time.
const INLINE_CSS: &str = r#"
/* === agogo.css (Sphinx agogo theme) === */
* { margin: 0px; padding: 0px; }
body { font-family: "Verdana", Arial, sans-serif; line-height: 1.4em; color: black; background-color: #eeeeec; }
a { color: #ce5c00; }
div.bodywrapper a, div.footer a { text-decoration: underline; }
.clearer { clear: both; }
.left { float: left; }
.right { float: right; }
img { border: 0; }

/* Header */
div.header, div.content, div.footer { width: 70em; margin-left: auto; margin-right: auto; }
div.header-wrapper { background: url(BGTOP_URI) top left repeat-x; border-bottom: 3px solid #2e3436; }
div.header { padding-top: 10px; padding-bottom: 10px; }
div.header .headertitle { font-family: "Georgia", "Times New Roman", serif; font-weight: normal; font-size: 180%; letter-spacing: .08em; margin-bottom: .8em; }
div.header .headertitle a { color: white; }
div.sidebar a, div.header a { text-decoration: none; }
div.sidebar a:hover, div.header a:hover { text-decoration: underline; }
p.logo { float: right; }
img.logo { border: 0; }

/* Headings */
h1, h2, h3, h4 { font-family: "Georgia", "Times New Roman", serif; font-weight: normal; color: #3465a4; margin-bottom: .8em; }
h1 { color: #204a87; }
h2 { padding-bottom: .5em; border-bottom: 1px solid #3465a4; }
a.headerlink { visibility: hidden; color: #dddddd; padding-left: .3em; }

/* Content */
div.content-wrapper { background-color: white; padding-top: 20px; padding-bottom: 20px; }
div.document { width: 50em; float: left; }
div.body { padding-right: 2em; text-align: justify; }
div.document h1 { line-height: 120%; }
div.document h2 { margin-top: .7em; }
div.document .section { margin-top: 1.7em; }
div.document .section:first-child { margin-top: 0px; }
img.align-center-qualimap { display: block; margin-left: auto; margin-right: auto; width: 100%; }

/* Sidebar */
div.sidebar { width: 20em; float: right; font-size: .9em; }
div.sidebar h3 { color: #2e3436; text-transform: uppercase; font-size: 130%; letter-spacing: .1em; }
div.sidebar ul { list-style-type: none; }
div.sidebar li.toctree-l1 a { display: block; padding: 1px; border: 1px solid #dddddd; background-color: #eeeeec; margin-bottom: .4em; padding-left: 3px; color: #2e3436; }
div.sidebar li.toctree-l1 a:hover { background-color: #dddddd; text-decoration: none; }
div.sidebar li.toctree-l1.current a { border-right: 5px solid #fcaf3e; font-weight: bold; }

/* Footer */
div.footer-wrapper { background: url(BGFOOTER_URI) top left repeat-x; border-top: 4px solid #babdb6; padding-top: 10px; padding-bottom: 10px; min-height: 80px; }
div.footer, div.footer a { color: #888a85; }
div.footer .right { text-align: right; }
div.footer .left { text-transform: uppercase; }

/* === report.css (Qualimap report additions) === */
table.summary { border: 0; width: 100%; vertical-align: top; }
div.table-summary { margin-left: auto; margin-right: auto; padding-bottom: 20px; }
td.column1 { width: 60%; }
table.hovertable { font-size: 14px; border: 1px solid #999999; border-collapse: collapse; }
table.hovertable th { background-color: #ffffff; border: 1px solid #a9c6c9; padding: 8px; }
table.hovertable tr { background-color: #ffffff; }
table.hovertable td { border: 1px solid #a9c6c9; padding: 8px; }
"#;

// ===================================================================
// HTML report data structure
// ===================================================================

/// All data needed to render the HTML report.
#[derive(Debug)]
pub struct ReportData<'a> {
    /// Sample/BAM stem name.
    pub sample_name: &'a str,
    /// Original BAM file path.
    pub bam_path: &'a str,
    /// GTF file path.
    pub gtf_path: &'a str,
    /// Whether the data is paired-end.
    pub paired: bool,
    /// Strandedness.
    pub stranded: Strandedness,
    /// Left-mapped proper pairs (PE only).
    pub left_proper: u64,
    /// Right-mapped proper pairs (PE only).
    pub right_proper: u64,
    /// Both-mapped proper pairs (PE only).
    pub both_proper: u64,
    /// Total read count (kept for API completeness).
    #[allow(dead_code)]
    pub read_count: u64,
    /// Primary alignments.
    pub primary_alignments: u64,
    /// Secondary alignments.
    pub secondary_alignments: u64,
    /// Non-unique alignments.
    pub alignment_not_unique: u64,
    /// Mapped reads on a sequence with no annotated exon.
    pub seq_not_found: u64,
    /// Reads aligned to genes (exonic).
    pub exonic_reads: u64,
    /// Ambiguous reads.
    pub ambiguous_reads: u64,
    /// No-feature reads.
    pub no_feature: u64,
    /// Not-aligned reads.
    pub not_aligned: u64,
    /// Intronic reads.
    pub intronic_reads: u64,
    /// Intergenic reads.
    pub intergenic_reads: u64,
    /// Reads overlapping exon boundary.
    pub overlapping_exon_reads: u64,
    /// Strand-specificity estimation: forward.
    pub ssp_fwd: u64,
    /// Strand-specificity estimation: reverse.
    pub ssp_rev: u64,
    /// Reads at splice junctions.
    pub reads_at_junctions: u64,
    /// Junction motif counts (e.g., "GT/AG" -> count).
    pub junction_motifs: &'a HashMap<String, u64>,
    /// 5' coverage bias.
    pub five_bias: f64,
    /// 3' coverage bias.
    pub three_bias: f64,
    /// 5'-3' coverage bias.
    pub five_three_bias: f64,
    /// Optional junction event counts: (known, partly_known, novel).
    pub junction_counts: Option<(u64, u64, u64)>,
}

// ===================================================================
// Public API
// ===================================================================

/// Write the self-contained `qualimapReport.html` to the output directory.
pub fn write_html_report(data: &ReportData, output_dir: &Path) -> Result<()> {
    let path = output_dir.join("qualimapReport.html");

    let html = render_report(data)?;

    let mut f = std::fs::File::create(&path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    f.write_all(html.as_bytes())
        .with_context(|| format!("Failed to write {}", path.display()))?;

    debug!("Wrote HTML report: {}", path.display());
    Ok(())
}

// ===================================================================
// Rendering
// ===================================================================

/// Render the full HTML report as a `String`.
fn render_report(data: &ReportData) -> Result<String> {
    let mut html = String::with_capacity(32 * 1024);

    // Build CSS with embedded data URIs for background images
    let css = INLINE_CSS
        .replace("BGTOP_URI", &format!("data:image/png;base64,{}", BGTOP_B64))
        .replace(
            "BGFOOTER_URI",
            &format!("data:image/png;base64,{}", BGFOOTER_B64),
        );

    // DOCTYPE + head
    html.push_str("<!DOCTYPE html>\n<html>\n<head>\n");
    html.push_str("<meta charset=\"utf-8\" />\n");
    html.push_str(&format!(
        "<title>Qualimap report: RNA Seq QC - {}</title>\n",
        data.sample_name
    ));
    html.push_str("<style>\n");
    html.push_str(&css);
    html.push_str("\n</style>\n");
    html.push_str("</head>\n<body>\n");

    // --- Header (matches Qualimap exactly) ---
    html.push_str("<div class=\"header-wrapper\">\n<div class=\"header\">\n");
    html.push_str(&format!(
        "<p class=\"logo\"><a href=\"http://qualimap.bioinfo.cipf.es/\">\
         <img class=\"logo\" src=\"data:image/png;base64,{}\" alt=\"Logo\"/></a></p>\n",
        LOGO_B64
    ));
    html.push_str(
        "<div class=\"headertitle\"><a href=\"\">Qualimap Report: RNA Seq QC</a></div>\n",
    );
    html.push_str("<div style=\"clear: both;\"></div>\n");
    html.push_str("</div>\n</div>\n");

    // --- Content wrapper ---
    html.push_str("<div class=\"content-wrapper\">\n<div class=\"content\">\n");
    html.push_str("<div class=\"document\">\n<div class=\"documentwrapper\">\n");
    html.push_str("<div class=\"bodywrapper\">\n<div class=\"body\">\n");

    // Section 1: Input data and parameters
    write_input_section(&mut html, data);

    // Section 2: Summary
    write_summary_section(&mut html, data);

    // Section 3+: Image sections
    write_image_sections(&mut html, data);

    // Close body/document
    html.push_str("</div>\n</div>\n</div>\n</div>\n");

    // --- Sidebar ---
    write_sidebar(&mut html, data);

    html.push_str("<div style=\"clear: both;\"></div>\n");
    html.push_str("</div>\n</div>\n");

    // --- Footer (matches Qualimap layout: date on left, generator on right) ---
    let date_str = format_local_datetime();
    html.push_str("<div class=\"footer-wrapper\">\n<div class=\"footer\">\n");
    html.push_str(&format!(
        "<div class=\"left\"><div class=\"footer\">{}</div></div>\n",
        date_str
    ));
    html.push_str(
        "<div class=\"right\"><div class=\"footer\">Generated by \
         <a href=\"https://github.com/seqeralabs/RustQC\">RustQC</a>, \
         compatible with QualiMap v.2.3</div></div>\n",
    );
    html.push_str("<div class=\"clearer\"></div>\n");
    html.push_str("</div>\n</div>\n");

    html.push_str("</body>\n</html>\n");
    Ok(html)
}

// -------------------------------------------------------------------
// Section writers
// -------------------------------------------------------------------

/// Section 1: Input data and parameters (matches Qualimap heading structure).
fn write_input_section(html: &mut String, data: &ReportData) {
    html.push_str("<div class=\"section\">\n");
    html.push_str(
        "<h2>Input data and parameters<a class=\"headerlink\" name=\"input\">&nbsp;</a></h2>\n",
    );
    html.push_str(
        "<p>Report generated by <a href=\"https://github.com/seqeralabs/RustQC\">RustQC</a>, \
         a Rust performance-rewrite of Qualimap. Results should be identical to Qualimap v2.3 \
         output.</p>\n",
    );

    // Qualimap wraps this in <div class="summary"><div class="table-summary"><h3>Input</h3>
    html.push_str("<div class=\"summary\">\n<div class=\"table-summary\">\n");
    html.push_str("<h3>Input</h3>\n");
    html.push_str("<table class=\"summary hovertable\">\n");

    let protocol = match data.stranded {
        Strandedness::Forward => "strand-specific-forward",
        Strandedness::Reverse => "strand-specific-reverse",
        Strandedness::Unstranded => "non-strand-specific",
    };

    // Date row — matches Qualimap's "Analysis date: " label
    let date_str = format_local_datetime();

    let rows = [
        ("Analysis date: ", date_str),
        ("BAM file: ", data.bam_path.to_string()),
        ("Counting algorithm:", "uniquely-mapped-reads".to_string()),
        ("GTF file: ", data.gtf_path.to_string()),
        (
            "Number of bases for 5'-3' bias computation:",
            format_long(NUM_PRIME_BASES as u64),
        ),
        (
            "Number of transcripts for 5'-3' bias computation:",
            format_long(NUM_TRANSCRIPTS_FOR_BIAS as u64),
        ),
        (
            "Paired-end sequencing:",
            if data.paired {
                "yes".to_string()
            } else {
                "no".to_string()
            },
        ),
        ("Protocol: ", protocol.to_string()),
        ("Sorting performed:", "no".to_string()),
    ];

    for (label, value) in &rows {
        table_row(html, label, value);
    }

    html.push_str("</table>\n</div>\n</div>\n");
    html.push_str("</div>\n");
}

/// Section 2: Summary (Reads alignment, Genomic origin, Coverage profile, Junction analysis).
fn write_summary_section(html: &mut String, data: &ReportData) {
    html.push_str("<div class=\"section\">\n");
    html.push_str("<h2>Summary<a class=\"headerlink\" name=\"summary\">&nbsp;</a></h2>\n");

    // --- Reads alignment (RNASeqQCAnalysis.prepareHtmlSummary) ---
    let mut alignment: Vec<(&str, String)> = Vec::new();
    if data.paired {
        alignment.push((
            "Number of mapped reads (left/right):",
            format!(
                "{} / {}",
                format_long(data.left_proper),
                format_long(data.right_proper)
            ),
        ));
        alignment.push((
            "Number of aligned pairs (without duplicates):",
            format_long(data.both_proper / 2),
        ));
    } else {
        alignment.push((
            "Number of mapped reads:",
            format_long(data.primary_alignments),
        ));
    }
    alignment.extend([
        (
            "Total number of alignments:",
            format_long(data.primary_alignments + data.secondary_alignments),
        ),
        (
            "Number of secondary alignments:",
            format_long(data.secondary_alignments),
        ),
        (
            "Number of non-unique alignments:",
            format_long(data.alignment_not_unique),
        ),
        ("Aligned to genes:", format_long(data.exonic_reads)),
        ("Ambiguous alignments:", format_long(data.ambiguous_reads)),
        ("No feature assigned:", format_long(data.no_feature)),
    ]);
    if data.seq_not_found > 0 {
        alignment.push((
            "Missing chromosome in annotation: ",
            format_long(data.seq_not_found),
        ));
    }
    alignment.push(("Not aligned:", format_long(data.not_aligned)));
    if data.stranded == Strandedness::Unstranded {
        let (fwd, rev) = super::strand_estimation(data.ssp_fwd, data.ssp_rev);
        alignment.push((
            "Strand specificity estimation (fwd/rev):",
            format!("{} / {}", format_decimal(fwd), format_decimal(rev)),
        ));
    }
    summary_table(html, "Reads alignment", &alignment);

    // --- Reads genomic origin ---
    // Percentages of exonic + "no feature" reads, as in the text report.
    let total_reads = data.exonic_reads + data.no_feature;
    let origin_row = |count: u64| {
        format!(
            "{} / {}",
            format_long(count),
            format_percentage(percentage_of(count, total_reads))
        )
    };
    summary_table(
        html,
        "Reads genomic origin",
        &[
            ("Exonic: ", origin_row(data.exonic_reads)),
            ("Intronic: ", origin_row(data.intronic_reads)),
            ("Intergenic: ", origin_row(data.intergenic_reads)),
            (
                "Intronic/intergenic overlapping exon: ",
                origin_row(data.overlapping_exon_reads),
            ),
        ],
    );

    // --- Transcript coverage profile ---
    summary_table(
        html,
        "Transcript coverage profile",
        &[
            ("5' bias:", format_decimal(data.five_bias)),
            ("3' bias:", format_decimal(data.three_bias)),
            ("5'-3' bias:", format_decimal(data.five_three_bias)),
        ],
    );

    // --- Junction analysis ---
    let mut junctions: Vec<(&str, String)> =
        vec![("Reads at junctions:", format_long(data.reads_at_junctions))];
    for (motif, pct) in super::junction_motif_rows(data.junction_motifs, data.reads_at_junctions) {
        junctions.push((motif, format_percentage(pct)));
    }
    summary_table(html, "Junction analysis", &junctions);

    html.push_str("</div>\n");
}

/// One `<h3>` summary table, as Qualimap's `HtmlReportGenerator.appendTableFromStats`.
fn summary_table(html: &mut String, title: &str, rows: &[(&str, String)]) {
    html.push_str("<div class=\"summary\">\n<div class=\"table-summary\">\n");
    html.push_str(&format!("<h3>{title}</h3>\n"));
    html.push_str("<table class=\"summary hovertable\">\n");
    for (label, value) in rows {
        table_row(html, label, value);
    }
    html.push_str("</table>\n</div>\n</div>\n");
}

/// Sections 3+: Image plots (matches Qualimap section structure exactly).
fn write_image_sections(html: &mut String, data: &ReportData) {
    // Each image gets its own <div class="section"> with h2 + headerlink anchor
    let image_sections = [
        ("Reads Genomic Origin", "Reads Genomic Origin"),
        (
            "Coverage Profile Along Genes (Total)",
            "Coverage Profile Along Genes (Total)",
        ),
        (
            "Coverage Profile Along Genes (Low)",
            "Coverage Profile Along Genes (Low)",
        ),
        (
            "Coverage Profile Along Genes (High)",
            "Coverage Profile Along Genes (High)",
        ),
        // Qualimap: h2 text is "Coverage Histogram (0-50X)" but anchor is "Transcript coverage histogram"
        (
            "Transcript coverage histogram",
            "Coverage Histogram (0-50X)",
        ),
    ];

    for (anchor, title) in &image_sections {
        // Determine the image filename — for the histogram, the file is named
        // "Transcript coverage histogram.png" to match Qualimap's output filename
        let img_name = if *anchor == "Transcript coverage histogram" {
            "Transcript coverage histogram"
        } else {
            title
        };
        html.push_str("<div class=\"section\">\n");
        html.push_str(&format!(
            "<h2>{}<a class=\"headerlink\" name=\"{}\">&nbsp;</a></h2>\n",
            title, anchor
        ));
        html.push_str(&format!(
            "<div><img width=\"100%\" src=\"images_qualimapReport/{}.png\"></div>\n",
            img_name
        ));
        html.push_str("</div>\n");
    }

    // Junction Analysis — only if non-zero junction counts exist
    if data.junction_counts.is_some_and(|(k, p, n)| k + p + n > 0) {
        html.push_str("<div class=\"section\">\n");
        html.push_str(
            "<h2>Junction Analysis<a class=\"headerlink\" name=\"Junction Analysis\">&nbsp;</a></h2>\n",
        );
        html.push_str(
            "<div><img width=\"100%\" src=\"images_qualimapReport/Junction Analysis.png\"></div>\n",
        );
        html.push_str("</div>\n");
    }
}

/// Sidebar with table of contents (matches Qualimap exactly — bare li items, no ul wrapper).
fn write_sidebar(html: &mut String, data: &ReportData) {
    html.push_str("<div class=\"sidebar\">\n");
    html.push_str("<h3>Contents</h3>\n");

    // Qualimap uses bare <li> items without a <ul> wrapper
    let toc_entries = [
        ("input", "Input data &amp; parameters"),
        ("summary", "Summary"),
        ("Reads Genomic Origin", "Reads Genomic Origin"),
        (
            "Coverage Profile Along Genes (Total)",
            "Coverage Profile Along Genes (Total)",
        ),
        (
            "Coverage Profile Along Genes (Low)",
            "Coverage Profile Along Genes (Low)",
        ),
        (
            "Coverage Profile Along Genes (High)",
            "Coverage Profile Along Genes (High)",
        ),
        (
            "Transcript coverage histogram",
            "Coverage Histogram (0-50X)",
        ),
    ];

    for (anchor, label) in &toc_entries {
        html.push_str(&format!(
            "<li class=\"toctree-l1\"><a href=\"#{}\">{}</a></li>\n",
            anchor, label
        ));
    }

    if data.junction_counts.is_some_and(|(k, p, n)| k + p + n > 0) {
        html.push_str(
            "<li class=\"toctree-l1\"><a href=\"#Junction Analysis\">Junction Analysis</a></li>\n",
        );
    }

    html.push_str("</div>\n");
}

// -------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------

/// Write a two-column table row with hover effect (matches Qualimap's exact colors).
fn table_row(html: &mut String, label: &str, value: &str) {
    html.push_str(
        "<tr onmouseover=\"this.style.backgroundColor='#EEEEEC';\" \
         onmouseout=\"this.style.backgroundColor='#FFFFFF';\">\n",
    );
    html.push_str(&format!(
        "<td class=\"column1\">{}</td><td class=\"column2\">{}</td>\n",
        label, value
    ));
    html.push_str("</tr>\n");
}

/// Format the current local date/time as `YYYY/MM/DD HH:MM:SS`.
///
/// Uses `time` crate's `OffsetDateTime` (already an indirect dependency via `env_logger`).
fn format_local_datetime() -> String {
    // Fallback: shell out to `date` command for reliable local time formatting
    // without adding a new dependency.
    if let Ok(output) = std::process::Command::new("date")
        .arg("+%Y/%m/%d %H:%M:%S")
        .output()
    {
        if output.status.success() {
            return String::from_utf8_lossy(&output.stdout).trim().to_string();
        }
    }
    // Ultimate fallback: UTC via std
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    // Simple UTC breakdown (no timezone offset)
    let days = secs / 86400;
    let time_of_day = secs % 86400;
    let h = time_of_day / 3600;
    let m = (time_of_day % 3600) / 60;
    let s = time_of_day % 60;
    // Compute date from days since epoch (1970-01-01)
    let (y, mo, d) = days_to_ymd(days);
    format!("{y:04}/{mo:02}/{d:02} {h:02}:{m:02}:{s:02}")
}

/// Convert days since Unix epoch to (year, month, day).
fn days_to_ymd(days: u64) -> (u64, u64, u64) {
    // Civil calendar algorithm from Howard Hinnant
    let z = days + 719468;
    let era = z / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    (y, mo, d)
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_table_row_html() {
        let mut html = String::new();
        table_row(&mut html, "Test Label", "Test Value");
        assert!(html.contains("Test Label"));
        assert!(html.contains("Test Value"));
        assert!(html.contains("onmouseover"));
        assert!(html.contains("#EEEEEC")); // Qualimap hover color
    }
}
