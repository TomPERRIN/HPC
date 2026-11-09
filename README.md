# CHP — TP3 : LS#3 3D sampling & Confidence Intervals

Benchmark de générateurs de nombres pseudo-aléatoires, puis Monte Carlo 3D
avec intervalles de confiance à 95 % (loi de Student).

## Fichiers

| Fichier | Rôle | Provenance |
|---|---|---|
| `OMP.c` | Exemple 2D fourni : estimation de PI par Monte Carlo, réplicats parallèles OpenMP, générateur xorshift64 thread-safe (état local par réplicat) | Énoncé du TP |
| `mt19937ar.c` / `mt19937ar.h` | Mersenne Twister MT19937 : `init_genrand(seed)` pour l'initialisation, `genrand_real2()` pour un tirage dans [0,1) | Code original de Matsumoto (Hiroshima), téléchargé |
| `MRG.h` | Paramètres (defines) des 3 générateurs de L'Ecuyer + prototypes | Générateurs de L'Ecuyer (code Moodle), restructuré : defines séparés dans le `.h` |
| `MRG.c` | MRG32k3a, MRG32k5a, MRG63k3a + **fonctions d'initialisation** (`MRGxx_init(seed)`, l'état tout-zéro est interdit) | Code Moodle, restructuré : un seul jeu de defines, états `static` renommés par générateur, init ajoutée |
| `benchmark.c` | Benchmark : N appels (défaut 10⁹) de chaque générateur dans [0,1), mesure du temps, ns/call et Mcalls/s ; une somme "sink" empêche l'optimiseur de supprimer la boucle | Écrit pour ce TP |
| `mt19937ar_ctx.c` / `.h` | Version **thread-safe à contexte** de MT19937 : l'état (`mt[624]` + `mti`) est dans une struct passée en paramètre au lieu de variables globales → chaque réplicat OpenMP possède son flux indépendant. Algorithme strictement inchangé (licence BSD conservée) | Dérivé du code original de Matsumoto |
| `mrg32k3a_ctx.c` / `.h` | Version **thread-safe à contexte** de MRG32k3a : les 6 doubles d'état dans une struct, paramètres réutilisés de `MRG.h`, algorithme inchangé | Dérivé du code de L'Ecuyer (Moodle) |
| `OMP3D.c` | **Partie 2** : Monte Carlo 3D — volume d'une sphère (4/3 π R³, R=1) ; n réplicats saisis par l'utilisateur, 10⁹ points/réplicat par défaut, choix du générateur (MT ou MRG32k3a) ; moyenne X̄(n), variance corrigée S²(n)/(n−1), erreur standard S(n)/√n, rayon de confiance R = t(n−1, 0.975) × S(n)/√n (table de Student intégrée), CI 95 % | Écrit pour ce TP |
| `LS # 3 - 3D sampling and Conf Intervals.pdf` | Sujet du TP | Énoncé |
| `Makefile` | `make` : tout compiler en -O2 (mcPI, bench_O2, mc3d) · `make bench` : + bench_O0 (sans optimisation) · `make clean` : nettoyer | Écrit pour ce TP |

## Compilation et exécution

```bash
make              # tout compiler en -O2 : mcPI + bench_O2
make bench        # + bench_O0 (sans optimisation) pour la comparaison
make clean        # supprimer .o et executables
```

Exemple 2D (fourni) :

```bash
make
OMP_NUM_THREADS=4 ./mcPI
```

Monte Carlo 3D + intervalles de confiance (partie 2) :

```bash
make                # compile mc3d
./mc3d              # saisie interactive : n réplicats, points/réplicat (0 = 10^9), générateur
# ex. 10 réplicats, 10^9 points, MT19937 :   10 / 0 / 1
# ex. 40 réplicats, 10^9 points, MRG32k3a :  40 / 0 / 2
```

Benchmark des générateurs :

```bash
./bench_O2                 # 10^9 appels par générateur (défaut), ~1 min 30
./bench_O2 10000000        # version rapide (10^7 appels), ~1 s
./bench_O0                 # même chose sans optimisation
```

## Benchmark des générateurs (10⁹ appels, en [0,1))

| Générateur | sans `-O` | avec `-O2` | gain `-O2` | débit `-O2` |
|---|---|---|---|---|
| `rand()` (référence) | 12.96 s (12.96 ns/call) | 12.91 s (12.91 ns/call) | ×1.0 | 77 Mcalls/s |
| **MT19937** (Matsumoto) | 4.41 s (4.41 ns/call) | **2.55 s (2.55 ns/call)** | **×1.7** | **392 Mcalls/s** |
| **MRG32k3a** (L'Ecuyer) | 26.74 s (26.74 ns/call) | 20.40 s (20.40 ns/call) | ×1.3 | 49 Mcalls/s |
| MRG32k5a (bonus) | 37.77 s (37.77 ns/call) | 28.72 s (28.72 ns/call) | ×1.3 | 35 Mcalls/s |
| MRG63k3a (bonus) | 24.86 s (24.86 ns/call) | 20.25 s (20.25 ns/call) | ×1.2 | 49 Mcalls/s |

### Observations

- **`rand()` est lent ET statistiquement faible** : ~13 ns/call, et `-O2` ne
  change rien (fonction dans la libc, déjà compilée). À bannir du calcul
  scientifique.
- **MT19937 est le plus rapide** des générateurs de qualité : ~5× plus rapide
  que MRG32k3a. Le gain de `-O2` vient de l'inlining de `genrand_int32`
  (mentionné sur la page de Matsumoto).
- **MRG32k3a est ~5× plus lent** : deux composantes récursives et un modulo
  flottant par tirage, mais il offre des flux multiples bien étudiés.
- MRG32k5a et MRG63k3a : haute qualité mais plus lents, comme annoncé dans
  le sujet.
- **À retenir** : toujours compiler avec `-O2` pour la simulation
  scientifique ; ne pas utiliser les `.o` d'une compilation optimisée pour
  produire la version sans optimisation.

## Monte Carlo 3D & intervalles de confiance (partie 2)

- `OMP3D.c` : Monte Carlo 3D — estimation du volume d'une sphère
  (4/3 π R³), 10⁹ points par réplicat, nombre de réplicats saisi par
  l'utilisateur, flux indépendant par réplicat.
- Intervalles de confiance à 95 % : moyenne X̄(n), variance corrigée
  S²(n)/(n−1), erreur standard S(n)/√n, rayon R = t(n−1, 0.975) × S(n)/√n
  (table de Student de l'annexe), comparaison 10 vs 40 réplicats, tests
  avec MT19937 et MRG32k3a.

### Résultats (10⁹ points par réplicat, 12 threads, sphère R = 1)

Volume exact : 4/3 π = 4.1887902048 · PI exact = 3.1415926536

| Configuration | Volume estimé | Erreur standard | t(n−1, 0.975) | Rayon de confiance R | Exact dans la CI ? | Temps |
|---|---|---|---|---|---|---|
| MT19937, n = 10 | 4.1887917688 | 3.38e-05 | 2.228 | 7.52e-05 | **oui** | 12.6 s |
| MT19937, n = 40 | 4.1887913002 | 2.13e-05 | 2.021 | 4.31e-05 | **oui** | 49.1 s |
| MRG32k3a, n = 10 | 4.1887618368 | 5.81e-05 | 2.228 | 1.29e-04 | **oui** | 68.1 s |
| MRG32k3a, n = 40 | 4.1887746644 | 2.16e-05 | 2.021 | 4.37e-05 | **oui** | 269.5 s |

Estimations de PI correspondantes (X̄ × 3/4) :

| Configuration | PI estimé | CI 95 % sur PI |
|---|---|---|
| MT19937, n = 10 | 3.1415938266 | [3.1415373943 , 3.1416502589] |
| MT19937, n = 40 | 3.1415934752 | [3.1415611209 , 3.1416258294] |
| MRG32k3a, n = 10 | 3.1415713776 | [3.1414742750 , 3.1416684802] |
| MRG32k3a, n = 40 | 3.1415809983 | [3.1415482533 , 3.1416137433] |

### Observations

- **La valeur exacte est dans l'intervalle de confiance à 95 % pour les
  4 expériences** : la méthode de l'annexe (Student) est validée.
- **Augmenter n de 10 à 40 réplicats réduit le rayon de confiance**
  (MT : 7.52e-05 → 4.31e-05 ; MRG : 1.29e-04 → 4.37e-05), pour deux
  raisons cumulées : t(n−1, 0.975) diminue (2.228 → 2.021) et l'erreur
  standard S(n)/√n diminue en ~1/√n.
- La variance estimée S²(n) fluctue elle-même d'une série de réplicats
  à l'autre (ex. MT : 1.14e-08 à n=10 contre 1.82e-08 à n=40) : le rayon
  de confiance ne suit donc pas exactement le facteur théorique √(10/40),
  mais l'ordre de grandeur reste conforme.
- **MT19937 ≈ 5× plus rapide que MRG32k3a** dans ce programme aussi
  (12.6 s contre 68.1 s pour 10 réplicats), cohérent avec le benchmark
  (2.5 ns/call contre 20.4 ns/call).
