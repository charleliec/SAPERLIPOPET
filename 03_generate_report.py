import os
import shutil
import tempfile
from datetime import datetime
from quarto import render

# Temporary isolated folder
date_time_launch = datetime.now().strftime("%d-%m-%Y_%Hh%Mmin%Ss")
work_dir = f"/tmp/cache_saperlipopet_generate_report_{date_time_launch}"
os.makedirs(work_dir, exist_ok=True)
os.chdir(work_dir)

# Isolate all caches and runtime
home = os.path.join(work_dir, "home")
os.makedirs(home, exist_ok=True)
os.environ["HOME"] = home
os.environ["XDG_CACHE_HOME"] = os.path.join(home, ".cache")
os.environ["DENO_DIR"] = os.path.join(home, ".deno_cache")
os.environ["QUARTO_USER_CACHE_DIR"] = os.path.join(home, ".quarto_cache")
os.environ["XDG_RUNTIME_DIR"] = os.path.join(home, ".runtime")

for d in ["XDG_CACHE_HOME", "DENO_DIR", "QUARTO_USER_CACHE_DIR", "XDG_RUNTIME_DIR"]:
    os.makedirs(os.environ[d], exist_ok=True)

# Copy inputs
qmd_file = "04_ref_eval_dashboard.qmd" if snakemake.params["reference_evaluation"] else "05_annotation_dashboard.qmd"

print("Generating the following report : ", qmd_file)

files_to_copy = [os.path.join(snakemake.params["scripts_dir"], qmd_file), os.path.join(snakemake.params["scripts_dir"], "Saperlipopet_logo.png")] + list(snakemake.input)
for f in files_to_copy:
    shutil.copy(f, work_dir)

# Parameters
params = {
        "REF_PATH": str(snakemake.input["ref"]),
        "QUERY_PATH": str(snakemake.input["query"]),
        "ANNOTATIONS_FILE": str(snakemake.input["annotations_file"]),
        "LOG_FILE_PATH": str(snakemake.log[0]),
    }

# Render
output_file = "dashboard.html"
render(
    input=qmd_file,
    output_file=output_file,
    #output_format="html",
    execute=True,
    execute_params=params,
    cache=False
)

# Move final report out
shutil.copy(output_file, snakemake.output[0])
shutil.copy(output_file, snakemake.output[0])

print("Report successfully generated:", snakemake.output[0])
