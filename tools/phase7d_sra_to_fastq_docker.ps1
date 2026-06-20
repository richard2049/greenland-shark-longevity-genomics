param(
    [string[]]$Runs = @("SRR32965275", "SRR32965277", "SRR32965276"),
    [int]$Threads = 8
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path ".").Path
$rawDir = Join-Path $repo "data\raw\rnaseq\SMIC_RETINA_PRJNA1246101_2026"
$logDir = Join-Path $repo "logs\rnaseq"
New-Item -ItemType Directory -Force $rawDir, $logDir | Out-Null

foreach ($run in $Runs) {
    $log = Join-Path $logDir "phase7d_${run}_sra_to_fastq.log"
    $container = "green_shark_phase7d_sra_$run"

    Write-Host "Starting $run; log: $log"
    $previousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & docker --context desktop-linux run --rm --name $container `
        -v "${repo}:/work" `
        -w /work `
        ncbi/sra-tools `
        sh /work/tools/phase7d_sra_to_fastq_container.sh $run $Threads *> $log
    $dockerExitCode = $LASTEXITCODE
    $ErrorActionPreference = $previousErrorActionPreference

    if ($dockerExitCode -ne 0) {
        Get-Content $log -Tail 120
        throw "SRA to FASTQ conversion failed for $run with exit code $dockerExitCode"
    }
    Get-Content $log -Tail 20
}
