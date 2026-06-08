# GitHub Pages Deployment Configuration for SynBio Pipeline

This directory contains the project website that can be deployed to GitHub Pages.

## Quick Deploy to GitHub Pages

### Option 1: Deploy from docs folder (Recommended)

```bash
# 1. Make sure you're on the main branch
git checkout main

# 2. Enable GitHub Pages in your repository settings:
#    - Go to Settings > Pages
#    - Source: Deploy from a branch
#    - Branch: main, folder: /docs
#    - Click Save

# 3. Your site will be live at: https://your-username.github.io/synbio-pipeline
```

### Option 2: Deploy to gh-pages branch

```bash
# Install gh-pages if you haven't already
npm install -g gh-pages

# Deploy the docs folder to gh-pages branch
gh-pages -d docs

# Your site will be live at: https://your-username.github.io/synbio-pipeline
```

### Option 3: Using GitHub Actions (Automatic)

The `.github/workflows/pages.yml` file automatically deploys to GitHub Pages when you push to main.

## File Structure

```
docs/
├── index.html              # Main HTML page
├── assets/
│   ├── css/
│   │   └── style.css      # All styles
│   ├── js/
│   │   └── main.js        # JavaScript functionality
│   └── images/            # (Optional) Add images here
└── README.md              # This file
```

## Customization

### Update Repository Links

In `index.html`, replace:
- `https://github.com/your-org/synbio-pipeline` with your actual repo URL
- `info@synbiopipeline.org` with your contact email

### Add Images

1. Place images in `docs/assets/images/`
2. Update references in `index.html`
3. For social media preview, create `og-image.png` (1200x630px)

### Modify Content

Edit `docs/index.html` to update:
- Hero text and subtitle
- Feature descriptions
- Organization information
- Footer links

## Local Development

You can preview the site locally:

```bash
# Using Python's built-in server
cd docs
python -m http.server 8000

# Visit http://localhost:8000
```

Or use any static file server:
- Live Server (VS Code extension)
- `npx serve docs`
- `php -S localhost:8000 -t docs`

## GitHub Pages Settings Checklist

- [ ] Repository is public (required for free GitHub Pages)
- [ ] GitHub Pages enabled in Settings > Pages
- [ ] Source set to "Deploy from a branch"
- [ ] Branch: main, Folder: /docs (or gh-pages branch)
- [ ] Custom domain configured (optional)
- [ ] HTTPS enforced (recommended)

## Troubleshooting

### Site not loading?
1. Check GitHub Actions tab for deployment errors
2. Verify file paths are relative (not absolute)
3. Clear browser cache

### CSS/JS not loading?
- Check browser console for 404 errors
- Verify paths in index.html are correct (relative to docs folder)

### Build fails?
- Ensure all files are committed
- Check workflow permissions in repository settings

## License

This website is part of the SynBio Pipeline project, released under the MIT License.
