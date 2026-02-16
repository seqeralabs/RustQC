// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";
import catppuccin from "@catppuccin/starlight";

// https://astro.build/config
export default defineConfig({
  site: process.env.SITE_URL || "https://ewels.github.io",
  base: process.env.BASE_PATH || "/RustQC",
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
          href: "https://github.com/ewels/RustQC",
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
          label: "Outputs",
          items: [
            { label: "dupRadar", slug: "outputs/dupradar" },
            { label: "featureCounts", slug: "outputs/featurecounts" },
            { label: "RSeQC", slug: "outputs/rseqc" },
            { label: "TIN & Gene Body", slug: "outputs/tin" },
            { label: "preseq", slug: "outputs/preseq" },
            { label: "Samtools", slug: "outputs/samtools" },
          ],
        },
        {
          label: "Benchmarks",
          items: [
            { label: "Combined", slug: "benchmarks/combined" },
            { label: "dupRadar", slug: "benchmarks/dupradar" },
            { label: "featureCounts", slug: "benchmarks/featurecounts" },
            { label: "RSeQC", slug: "benchmarks/rseqc" },
            { label: "preseq", slug: "benchmarks/preseq" },
            { label: "Samtools", slug: "benchmarks/samtools" },
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
