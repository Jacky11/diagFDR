// Browser smoke test of the Shiny app (Playwright). Milestone 5 check.
// 1. start the app:   Rscript -e 'DiffPTM::run_app(options = list(port = 3939))'
// 2. in a folder containing protein.tsv, ptm.tsv, annotation.tsv (Generic format):
//    npm i playwright && npx playwright install chromium && node e2e_smoke.js
// Produces shot*.png screenshots and app_results.xlsx; prints JS errors (expect "none").
const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const page = await b.newPage({ viewport: { width: 1400, height: 950 } });
  const errors = []; page.on('pageerror', e => errors.push(String(e)));
  await page.goto('http://127.0.0.1:3939/'); await page.waitForTimeout(2000);
  await page.evaluate(() => $('#import-software')[0].selectize.setValue('Generic'));
  await page.setInputFiles('#import-protein_file', 'protein.tsv');
  await page.setInputFiles('#import-ptm_file', 'ptm.tsv');
  await page.waitForTimeout(1500);
  await page.click('#import-go'); await page.waitForTimeout(5000);
  await page.screenshot({ path: 'shot1_import.png' });
  await page.click('text=2. Design'); await page.waitForTimeout(1000);
  await page.setInputFiles('#design-ann_file', 'annotation.tsv'); await page.waitForTimeout(2000);
  await page.click('#design-apply'); await page.waitForTimeout(3000);
  await page.screenshot({ path: 'shot2_design.png' });
  await page.click('text=3. Processing'); await page.waitForTimeout(1000);
  await page.check('#process-norm input[value="none"]'); await page.check('#process-impute input[value="none"]');
  await page.click('#process-run'); await page.waitForTimeout(15000);
  await page.screenshot({ path: 'shot3_process.png' });
  await page.click('text=4. Statistics'); await page.waitForTimeout(1000);
  await page.evaluate(() => $('#stats-pi0')[0].selectize.setValue('bh'));
  await page.click('#stats-run'); await page.waitForTimeout(20000);
  await page.click('text=5. Results'); await page.waitForTimeout(1500);
  await page.evaluate(() => $('#results-analysis')[0].selectize.setValue('ptm_DPU')); await page.waitForTimeout(8000);
  await page.screenshot({ path: 'shot5_results.png', fullPage: true });
  await page.click('text=6. Download'); await page.waitForTimeout(3000);
  await page.screenshot({ path: 'shot6_download.png' });
  const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 60000 }), page.click('#download-xlsx')]);
  await dl.saveAs('app_results.xlsx');
  console.log('JS errors:', errors.length ? errors : 'none');
  await b.close();
})().catch(e => { console.error('FAILED', e); process.exit(1); });
