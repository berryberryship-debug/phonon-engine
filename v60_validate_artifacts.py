import csv, os, sys
import numpy as np

CSV_PATH = "phonon_dispersion_501.csv"
PNG_PATH = "phonon_curves_501.png"
EXPECTED_HEADERS = ["q_m^-1", "omega_acoustic_rad_s^-1", "omega_optical_rad_s^-1", "vg_aco_ms", "vg_opt_ms"]

def fail(msg):
    print(f"[FAIL] {msg}")
    sys.exit(1)

def main():
    print("=== V60 : VALIDATION DES ARTEFACTS ===")
    
    if not os.path.exists(CSV_PATH): fail(f"Manquant : {CSV_PATH}")
    if not os.path.exists(PNG_PATH): fail(f"Manquant : {PNG_PATH}")
    
    png_size = os.path.getsize(PNG_PATH)
    if png_size < 10240: fail(f"PNG trop petit ({png_size} octets).")
    print(f"  [PASS] Fichier PNG valide ({png_size} octets)")

    with open(CSV_PATH, mode="r", newline="") as f:
        reader = csv.reader(f)
        headers = next(reader)
        rows = list(reader)
        
    if headers != EXPECTED_HEADERS: fail(f"En-tetes incorrects : {headers}")
    print("  [PASS] En-tetes CSV conformes")
    
    if len(rows) != 501: fail(f"Lignes : {len(rows)} (attendu: 501)")
    print("  [PASS] Nombre de lignes exact (501)")

    data = np.array(rows, dtype=float)
    if not np.all(np.isfinite(data)): fail("Valeurs NaN ou infinies.")
    print("  [PASS] Toutes les valeurs sont finies")

    if data[250, 0] != 0.0 or data[250, 1] != 0.0 or data[250, 2] <= 0.0:
        fail("Erreur physique au point q=0.")
    print("  [PASS] Coherence physique du point central")
    print(" SUCCES : Les artefacts sont valides.")

if __name__ == "__main__":
    main()
