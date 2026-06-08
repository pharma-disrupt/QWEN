# 🧬 SynBio Pipeline - GitHub Pages Project Site

## ✅ Files Created for Deployment

Your project website is ready to deploy! Here's what was created:

### Website Files
```
docs/
├── index.html              # Complete single-page website (20KB)
├── README.md               # Documentation for the site
├── assets/
│   ├── css/
│   │   └── style.css      # Modern responsive styles (18KB)
│   ├── js/
│   │   └── main.js        # Interactive functionality (3KB)
│   └── images/
│       └── README.md      # Image asset guidelines
│       └── .gitkeep       # Placeholder for future images
```

### GitHub Configuration
```
.github/workflows/
└── pages.yml              # Auto-deploy workflow
```

### Documentation
```
DEPLOYMENT.md              # Complete deployment guide
```

---

## 🚀 Deploy in 3 Steps

### Option 1: GitHub Pages from docs folder (Recommended)

1. **Go to your GitHub repo** → Settings → Pages
2. **Select:** Source = "Deploy from a branch", Branch = main, Folder = /docs
3. **Click Save** and push your code:

```bash
git add docs/ .github/ DEPLOYMENT.md
git commit -m "Add project website for GitHub Pages"
git push origin main
```

**Live URL:** `https://YOUR-USERNAME.github.io/YOUR-REPO`

### Option 2: Automatic via GitHub Actions

The `pages.yml` workflow automatically deploys when you push to main.

Just enable Actions in your repo settings, then:
```bash
git push origin main
```

### Option 3: Deploy to gh-pages branch

```bash
npm install -g gh-pages
gh-pages -d docs
```

---

## 🎨 Website Features

✅ **Modern Design** - Clean, professional UI with gradient accents  
✅ **Fully Responsive** - Works on desktop, tablet, and mobile  
✅ **5 Stage Overview** - Explains complete DBTL workflow  
✅ **Organism Cards** - Shows all 5 supported microorganisms  
✅ **Interactive Tabs** - Docker, pip, and source installation options  
✅ **Smooth Animations** - Scroll-triggered animations  
✅ **SEO Optimized** - Meta tags, Open Graph support  
✅ **Dark Navbar** - Glassmorphism effect on scroll  
✅ **No Dependencies** - Pure HTML/CSS/JS, no frameworks  

---

## 📱 Sections Included

1. **Hero** - Tagline, CTA buttons, key statistics
2. **Features** - 5 stages of the pipeline
3. **Organisms** - E. coli, Yeast, B. subtilis, C. glutamicum, P. putida
4. **How It Works** - Step-by-step workflow
5. **Installation** - Docker, pip, source code tabs
6. **Documentation** - Links to guides and API reference
7. **Footer** - Contact, social links, resources

---

## 🎯 Customization Tips

### Update Repository Links
Edit `docs/index.html`:
- Replace `https://github.com/your-org/synbio-pipeline` 
- Update email: `info@synbiopipeline.org`

### Add Custom Domain
Create `docs/CNAME`:
```
yourdomain.com
```

Then configure DNS with your domain provider.

### Add Images (Optional)
Place in `docs/assets/images/`:
- `og-image.png` (1200x630) - Social media preview
- Screenshots, logos, etc.

*Note: Site works perfectly with emoji icons only!*

---

## 🔍 Test Locally

Before deploying, preview locally:

```bash
cd docs
python -m http.server 8000
```

Visit: http://localhost:8000

---

## 📊 Post-Deployment Checklist

- [ ] Site loads at GitHub Pages URL
- [ ] All CSS styles apply correctly
- [ ] JavaScript tabs work
- [ ] Mobile view is responsive
- [ ] No console errors (F12)
- [ ] Smooth scrolling works
- [ ] All sections visible

---

## 🛠️ Troubleshooting

**404 Error?**
- Wait 1-2 minutes after pushing
- Verify GitHub Pages is enabled
- Check branch/folder settings

**Styles not loading?**
- Ensure paths are relative (not absolute)
- Check browser console for 404s

**Need help?**
See `DEPLOYMENT.md` for detailed troubleshooting guide.

---

## 📄 License

This website is part of SynBio Pipeline, released under MIT License.

---

**Ready to deploy? Run these commands:**

```bash
git add docs/ .github/ DEPLOYMENT.md
git commit -m "🚀 Add project website for GitHub Pages"
git push origin main
```

Then enable GitHub Pages in Settings → Pages → Deploy from branch: main, /docs

**Your site will be live in 2-3 minutes! 🎉**
