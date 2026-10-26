from ase.build import bulk

# Création d'un vrai cristal de Silicium (structure diamant)
# 5.43 est le paramètre de maille réel du Silicium en Angströms
silicium = bulk('Si', 'diamond', a=5.43, cubic=True)

print("=== CRISTAL 3D DE SILICIUM ===")
print(f"Nombre d'atomes : {len(silicium)}")
print(f"Formule chimique : {silicium.get_formula()}")

# Affichage des positions des 3 premiers atomes
print("\nPositions (x, y, z) des 3 premiers atomes :")
for i in range(3):
    print(f"Atome {i} : {silicium.positions[i]}")

# On sauvegarde le cristal dans un fichier standard pour la 3D
silicium.write('cristal_silicium.xyz')
print("\n[+] Cristal sauvegardé dans 'cristal_silicium.xyz'")
