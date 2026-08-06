# =====================================================================
#  build.ps1 — сборка учебного курса (зона D)
#
#  .\build.ps1            только LaTeX (быстро)
#  .\build.ps1 -Figures   сначала пересобрать все графики, потом LaTeX
#  .\build.ps1 -Clean     убрать промежуточные файлы сборки
#
#  Требуется MiKTeX (pdflatex, latexmk) и .venv проекта для графиков.
# =====================================================================
param(
    [switch]$Figures,
    [switch]$Clean
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$python = Join-Path $root "..\..\..\.venv\Scripts\python.exe"

if ($Clean) {
    Get-ChildItem $root -Include *.aux, *.log, *.out, *.toc, *.fls, *.fdb_latexmk -Recurse |
        Remove-Item -Force
    Write-Host "[build] промежуточные файлы удалены"
    exit 0
}

if ($Figures) {
    Write-Host "[build] пересборка графиков..."
    Push-Location (Join-Path $root "figures")
    try {
        Get-ChildItem -Filter "fig*.py" | ForEach-Object {
            & $python $_.Name
            if ($LASTEXITCODE -ne 0) { throw "график $($_.Name) не собрался" }
        }
    } finally { Pop-Location }
}

# Если main.pdf открыт в просмотрщике, Windows держит файл и pdflatex падает с
# «I can't write on file `main.pdf'». Закрывать чужую программу нельзя, поэтому
# в таком случае собираем в подкаталог build\ и сообщаем, где лежит результат.
$outDir = $root
$mainPdf = Join-Path $root "main.pdf"
if (Test-Path $mainPdf) {
    try {
        $fs = [System.IO.File]::Open($mainPdf, 'Open', 'Write', 'None')
        $fs.Close()
    } catch {
        $outDir = Join-Path $root "build"
        if (-not (Test-Path $outDir)) { New-Item -ItemType Directory $outDir | Out-Null }
        Write-Host "[build] main.pdf занят другой программой -> собираю в build\"
        # pdflatex ищет входные файлы сначала в текущем каталоге: оставшиеся в корне
        # main.toc/main.aux от прежней сборки перебили бы свежие из build\, и
        # оглавление собралось бы по устаревшим данным.
        Get-ChildItem $root -Include main.aux, main.toc, main.out, main.log -File |
            Remove-Item -Force -ErrorAction SilentlyContinue
    }
}

Write-Host "[build] компиляция main.tex..."
Push-Location $root
try {
    # Три прохода pdflatex, а НЕ latexmk: latexmk в этой установке MiKTeX
    # требует perl, которого в системе нет. Три прохода нужны, чтобы
    # устоялись оглавление и перекрёстные ссылки.
    for ($i = 1; $i -le 3; $i++) {
        & pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$outDir" main.tex | Out-Null
        if ($LASTEXITCODE -ne 0) {
            Write-Host "[build] ОШИБКА на проходе $i — выдержка из main.log:"
            Select-String -Path (Join-Path $outDir "main.log") -Pattern '^!' -Context 0, 5 |
                Select-Object -First 3 | ForEach-Object { $_.Line; $_.Context.PostContext }
            throw "LaTeX завершился с ошибкой (полный лог: $(Join-Path $outDir 'main.log'))"
        }
    }
} finally { Pop-Location }

$pdf = Join-Path $outDir "main.pdf"
$size = [math]::Round((Get-Item $pdf).Length / 1KB)
Write-Host "[build] готово: $pdf ($size КБ)"
