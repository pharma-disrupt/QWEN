# Deploy SynBio Pipeline Website to GitHub Pages

Complete guide to deploy your project website.

## 🚀 Quick Start (3 Steps)

### Step 1: Verify Repository Settings

```bash
# Make sure your repository is public
# GitHub Pages free tier requires public repositories

# Check your current branch
git branch

# Should be on main or master
git checkout main
```

### Step 2: Enable GitHub Pages

1. Go to your GitHub repository
2. Click **Settings** → **Pages** (left sidebar)
3. Under **Source**, select:
   - **Deploy from a branch**
   - Branch: `main` 
   - Folder: `/docs`
4. Click **Save**

### Step 3: Push Your Changes

```bash
# Add all new files
git add docs/ .github/

# Commit
git commit -m "Add project website for GitHub Pages"

# Push to main branch
git push origin main
```

Your site will be live at: `https://YOUR-USERNAME.github.io/YOUR-REPO`

---

## 📋 Pre-Deployment Checklist

- [ ] Repository is **public** (required for free GitHub Pages)
- [ ] All website files are in the `docs/` folder
- [ ] `docs/index.html` exists and is valid HTML
- [ ] Paths in HTML are relative (not absolute)
- [ ] GitHub Actions enabled in repository settings
- [ ] No sensitive data in website files

---

## 🔧 Alternative Deployment Methods

### Method A: Using GitHub Actions (Automatic)

The `.github/workflows/pages.yml` file automatically deploys when you push.

**Requirements:**
- GitHub Actions must be enabled
- Repository must allow Actions

**Enable Actions:**
1. Go to Settings → Actions
2. Select "Allow all actions"
3. Save changes

**Deploy:**
```bash
git push origin main
```

The workflow will automatically deploy to GitHub Pages.

### Method B: Deploy to gh-pages Branch

Use this if you want to keep website separate from code.

```bash
# Install gh-pages tool
npm install -g gh-pages

# Navigate to project root
cd /path/to/synbio-pipeline

# Deploy docs folder to gh-pages branch
gh-pages -d docs

# Wait ~1 minute, then visit:
# https://YOUR-USERNAME.github.io/YOUR-REPO
```

### Method C: Manual gh-pages Branch

```bash
# Create orphan gh-pages branch
git checkout --orphan gh-pages

# Remove all files except docs content
git rm -rf .
git checkout main -- docs/*
mv docs/* .
rm -rf docs

# Add and commit
git add .
git commit -m "Deploy website to gh-pages"

# Push gh-pages branch
git push origin gh-pages

# In GitHub Settings → Pages, select gh-pages branch
git checkout main
```

---

## ⚙️ Custom Domain Setup (Optional)

### Using a Custom Domain

1. Buy a domain (e.g., from Namecheap, GoDaddy)
2. In GitHub Settings → Pages → Custom domain:
   - Enter your domain (e.g., `synbiopipeline.org`)
   - Check "Enforce HTTPS"
3. Configure DNS records with your domain provider:

**For apex domain (example.org):**
```
Type: A
Name: @
Value: 185.199.108.153
TTL: 3600
```
(Create 4 A records with GitHub's IPs: 185.199.108.153, 185.199.109.153, 185.199.110.153, 185.199.111.153)

**For subdomain (www.example.org):**
```
Type: CNAME
Name: www
Value: YOUR-USERNAME.github.io
TTL: 3600
```

4. Create `docs/CNAME` file:
```
synbiopipeline.org
```

5. Commit and push:
```bash
echo "synbiopipeline.org" > docs/CNAME
git add docs/CNAME
git commit -m "Add custom domain"
git push
```

---

## 🐛 Troubleshooting

### Site Returns 404

**Causes:**
- GitHub Pages not enabled
- Wrong branch/folder selected
- Files not in correct location

**Solutions:**
1. Verify Settings → Pages shows "Your site is live"
2. Check branch is `main` and folder is `/docs` (or `gh-pages`)
3. Ensure `docs/index.html` exists
4. Wait 1-2 minutes after pushing

### CSS/JS Not Loading

**Check:**
1. Open browser DevTools (F12)
2. Look for 404 errors in Console
3. Verify paths are relative in `index.html`

**Fix paths if needed:**
```html
<!-- ❌ Wrong (absolute) -->
<link rel="stylesheet" href="/assets/css/style.css">

<!-- ✅ Correct (relative) -->
<link rel="stylesheet" href="assets/css/style.css">
```

### Build Fails in GitHub Actions

**Check:**
1. Go to Actions tab
2. Click failed workflow
3. Read error logs

**Common fixes:**
- Enable Actions in repository settings
- Check workflow permissions
- Verify `.github/workflows/pages.yml` syntax

### Changes Not Appearing

**Solutions:**
1. Hard refresh browser (Ctrl+Shift+R or Cmd+Shift+R)
2. Clear browser cache
3. Wait 2-5 minutes for GitHub to rebuild
4. Check Actions tab for deployment status

---

## 📊 Monitoring Deployment

### Check Deployment Status

1. **GitHub Actions Tab:** Shows build/deploy progress
2. **Settings → Pages:** Shows live URL and status
3. **Environment tab:** Shows deployment history

### View Deployment Logs

```bash
# In GitHub UI:
# Actions → Deploy to GitHub Pages → Latest run
```

---

## 🔒 Security Best Practices

- Never commit API keys or secrets to the website
- Use environment variables for sensitive config
- Keep dependencies updated
- Enable Dependabot alerts
- Review third-party scripts carefully

---

## 📈 Post-Deployment

### Share Your Site

- Add URL to repository description
- Update README with live demo link
- Share on social media
- Add to project documentation

### Analytics (Optional)

Add Google Analytics or similar:

```html
<!-- Add before </head> in docs/index.html -->
<script async src="https://www.googletagmanager.com/gtag/js?id=GA_MEASUREMENT_ID"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'GA_MEASUREMENT_ID');
</script>
```

---

## 📞 Need Help?

- **GitHub Pages Docs:** https://docs.github.com/en/pages
- **Community Forum:** https://github.community
- **This Project Issues:** https://github.com/YOUR-REPO/issues

---

## ✅ Success Indicators

You'll know it worked when:
- ✅ Settings → Pages shows green checkmark
- ✅ Site loads at `https://USERNAME.github.io/REPO`
- ✅ All styles and scripts load correctly
- ✅ Mobile responsive design works
- ✅ No console errors in browser

**Congratulations! Your SynBio Pipeline website is live! 🎉**
