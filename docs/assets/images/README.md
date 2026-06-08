# SynBio Pipeline - Project Website Assets

This directory contains placeholder information for project assets.

## Required Images

Create these images for a complete website:

### Open Graph / Social Media Preview
- **File:** `og-image.png`
- **Size:** 1200 x 630 pixels
- **Content:** SynBio Pipeline logo with tagline
- **Usage:** Social media sharing preview

### Logo (Optional)
- **File:** `logo.svg` or `logo.png`
- **Size:** Scalable (SVG recommended)
- **Content:** Project logo
- **Usage:** Navigation bar

### Screenshots (Optional)
- **File:** `screenshot-pipeline.png`
- **Size:** 1920 x 1080 pixels
- **Content:** Pipeline running in terminal
- **Usage:** Hero section or features

- **File:** `screenshot-results.png`
- **Size:** 1920 x 1080 pixels
- **Content:** Example output report
- **Usage:** Documentation section

## Creating Assets

### Quick Logo Creation
You can create a simple logo using:
```svg
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">
  <text y=".9em" font-size="90">🧬</text>
</svg>
```

### Using the DNA Emoji
The current site uses emoji (🧬) which works well and requires no assets.

## Asset Optimization

Before committing images:
1. Compress PNGs with tools like TinyPNG
2. Use WebP format for better compression
3. Keep total image size under 500KB per image
4. Use SVG for logos and icons when possible

## Directory Structure

```
docs/assets/
├── images/
│   ├── og-image.png      # Social media preview (create this)
│   ├── logo.svg          # Logo (optional)
│   └── screenshots/      # Optional screenshots
├── css/
│   └── style.css
└── js/
    └── main.js
```

## Current Status

✅ HTML structure complete
✅ CSS styling complete  
✅ JavaScript functionality complete
⏳ Images optional (emoji fallback works)

The website is fully functional without custom images - emojis are used as placeholders.
