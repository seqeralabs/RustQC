// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";
import catppuccin from "@catppuccin/starlight";

// https://astro.build/config
export default defineConfig({
  site: process.env.SITE_URL || "https://rustqc.netlify.app",
  integrations: [
    starlight({
      expressiveCode: {
        defaultProps: {
          frame: "none",
        },
      },
      title: "RustQC",
      logo: {
        src: "./src/assets/RustQC-icon.svg",
      },
      favicon: "/favicon.svg",
      social: [
        {
          icon: "github",
          label: "GitHub",
          href: "https://github.com/seqeralabs/RustQC",
        },
      ],
      sidebar: [
        {
          label: "Getting Started",
          items: [
            { label: "Introduction", slug: "getting-started/introduction" },
            { label: "Installation", slug: "getting-started/installation" },
            { label: "Quick Start", slug: "getting-started/quickstart" },
          ],
        },
        {
          label: "Usage",
          items: [
            { label: "CLI Reference", slug: "usage/cli-reference" },
            {
              label: "Configuration File",
              slug: "usage/configuration",
            },
          ],
        },
        {
          label: "RNA",
          items: [
            { label: "Benchmark Details", slug: "rna/benchmark-details" },
            { label: "dupRadar", slug: "rna/dupradar" },
            { label: "featureCounts", slug: "rna/featurecounts" },
            { label: "RSeQC", slug: "rna/rseqc" },
            { label: "Qualimap", slug: "rna/qualimap" },
            { label: "Preseq", slug: "rna/preseq" },
            { label: "Samtools", slug: "rna/samtools" },
            { label: "Performance & Tuning", slug: "rna/performance" },
          ],
        },
        {
          label: "About",
          items: [
            {
              label: "Credits & Citation",
              slug: "about/credits",
            },
            {
              label: "Contributing",
              slug: "about/contributing",
            },
          ],
        },
      ],
      customCss: ["./src/styles/custom.css"],
      plugins: [
        catppuccin({
          dark: { flavor: "mocha", accent: "peach" },
          light: { flavor: "latte", accent: "peach" },
        }),
      ],
    }),
  ],
});
