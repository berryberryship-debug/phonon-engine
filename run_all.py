import subprocess, sys, time

def run_script(script_name):
    print(f"\n{'-'*60}\nExecution de : {script_name}\n{'-'*60}")
    start = time.time()
    result = subprocess.run([sys.executable, script_name], capture_output=False)
    if result.returncode != 0:
        print(f"\nERREUR : {script_name} a echoue (code {result.returncode})")
        sys.exit(result.returncode)
    print(f"OK : {script_name} termine en {time.time()-start:.2f}s")

if __name__ == "__main__":
    print("="*60 + "\nLANCEMENT DU PIPELINE COMPLET\n" + "="*60)
    run_script("v57_production_501.py")
    run_script("v60_validate_artifacts.py")
    run_script("v58_regression_test.py")
    print("\n" + "="*60 + "\nSUCCES COMPLET : Generation, Artefacts et Physique valides.\n" + "="*60)
