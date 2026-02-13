// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";
import catppuccin from "@catppuccin/starlight";

// https://astro.build/config
export default defineConfig({
  site: "https://ewels.github.io",
  base: "/RustQC",
  integrations: [
    starlight({
      title: "RustQC",
      logo: {
        light: "./src/assets/RustQC-logo.svg",
        dark: "./src/assets/RustQC-logo-darkbg.svg",
        replacesTitle: true,
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
            { label: "Output Files", slug: "usage/output-files" },
          ],
        },
        {
          label: "Guide",
          items: [
            {
              label: "RNA-seq Duplicate Analysis",
              slug: "guide/rna-duprate",
            },
            {
              label: "Interpreting Plots",
              slug: "guide/interpreting-plots",
            },
            {
              label: "featureCounts Output",
              slug: "guide/featurecounts",
            },
          ],
        },
        {
          label: "Benchmarks",
          items: [
            {
              label: "Performance",
              slug: "benchmarks/performance",
            },
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
      plugins: [
        catppuccin({
          dark: { flavor: "mocha", accent: "peach" },
          light: { flavor: "latte", accent: "peach" },
        }),
      ],
    }),
  ],
});
