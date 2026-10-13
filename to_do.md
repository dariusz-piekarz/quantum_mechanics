# Plan rozwoju: orbitale cząsteczek i solver Kohna–Shama

## Cel i przyjęty zakres

Celem jest obliczanie stanów elektronowych **izolowanych cząsteczek** na trójwymiarowej siatce w przestrzeni rzeczywistej, a następnie ich analiza i wizualizacja. Pierwszy kompletny wynik powinien obejmować energię całkowitą, gęstość elektronową, orbitale i raport zbieżności. Rozwijamy najpierw układy jednoelektronowe, potem zamkniętopowłokowy KS-LDA dla He i H₂. Spin, PBE i większe cząsteczki przychodzą później.

Stan wyjściowy: istnieją `GridMatrix` z różnymi krokami osi, operator różnic skończonych, `eigen_solver`, potencjały Coulomba i harmoniczny oraz wizualizacja GLMakie. Są testy zbieżności energii swobodnej cząstki przy zagęszczaniu siatki. Nie ma jeszcze jawnego modelu cząsteczki, normalizacji orbitalu z objętością komórki, solvera Poissona, energii całkowitej KS ani pętli SCF. Aktualne izopowierzchnie używają indeksów tablicy jako współrzędnych.

Wszystkie poniższe wzory zakładają **jednostki atomowe**: odległość w bohr, energia w hartree, ładunek elektronu równy `-1`. Konwersje z Å i eV mają być jawne na wejściu/wyjściu. Nie mieszamy tych jednostek wewnątrz solvera.

## Granica między obliczeniami a wizualizacją

| Warstwa | Odpowiedzialność | Zależności |
| --- | --- | --- |
| `QuantumMechanics` — rdzeń | Geometria, siatka, potencjały, operatory, stany własne, gęstość, Poisson, SCF, energie, dane wynikowe i ich zapis do formatu wolumetrycznego. | Biblioteki numeryczne; **bez GLMakie, Meshing, GeometryBasics i LaTeXStrings**. |
| Opcjonalne rozszerzenie `QuantumMechanicsMakieExt` | Wykresy, izopowierzchnie, kolory fazy, osie, atomy, eksport PNG/SVG. | GLMakie, Meshing, GeometryBasics, LaTeXStrings. |
| `example/` | Skrypty łączące obie warstwy. | Przykłady numeryczne uruchamiają się bez środowiska graficznego. |

Preferowany mechanizm: rozszerzenie pakietu Julia (`ext/`, `[weakdeps]` i `[extensions]` w `Project.toml`). Jeśli kod graficzny urośnie znacznie bardziej niż rdzeń, można go później wydzielić do osobnego pakietu. Przy migracji zachować czytelne komunikaty dla dotychczasowych funkcji `plot_*` i zaktualizować przykłady.

**Kontrakt danych:** solver zwraca obiekt zawierający osie `x, y, z`, wartości pola w układzie `(nx, ny, nz)`, jednostkę, opis układu i metadane wyniku. Wizualizacja dostaje gotowe pole lub `SCFResult`; nie buduje hamiltonianu, nie uruchamia SCF i nie zmienia energii czy gęstości. Zapis pliku z próbkami pola, np. Gaussian cube, powinien być możliwy bez GLMakie.

## Kolejność prac

### 0. Oddzielić rdzeń od grafiki

- [x] Usunąć `using GLMakie`, `Meshing`, `GeometryBasics` i `LaTeXStrings` z `src/QuantumMechanics.jl`; przenieść kod `src/visualisation/` do opcjonalnego rozszerzenia.
- [x] Przenieść te pakiety z obowiązkowych zależności do zależności opcjonalnych. Zachować możliwość uruchamiania testów rdzenia bez OpenGL i bez wyświetlacza.
- [x] Rozdzielić przykłady na obliczeniowe i graficzne; opisać, co trzeba doinstalować do wykresów.
- [ ] Zapisać prosty kontrakt typów: `Grid3D`, `OrbitalSet`/`VolumetricField` i docelowo `SCFResult`. Rdzeń nie powinien zwracać obiektów Makie.

**Gotowe, gdy:** `using QuantumMechanics` oraz testy solvera działają bez instalacji GLMakie, a osobny test graficzny nadal tworzy dotychczasowe wykresy.

### 1. Uporządkować siatkę, geometrię i jednostki

- [ ] Wprowadzić `Grid3D` przechowujące jawnie trzy osie, `(nx, ny, nz)`, `(hx, hy, hz)`, początek układu i objętość komórki. Obecne `GridMatrix` może pozostać zgodnym interfejsem; nie trzeba trzymać całej macierzy `N × 3`, jeśli współrzędne można generować z osi.
- [ ] Ustalić jeden porządek spłaszczania: `x` najszybciej, następnie `y`, `z`. Dodać bezpieczne `reshape`/`flatten` i testy odwrotności.
- [ ] Dodać typ `Molecule`: liczby atomowe `Z_A`, pozycje jąder `R_A`, ładunek całkowity, liczba elektronów i krotność spinowa. Sprawdzać spójność tych danych.
- [ ] Jawnie oddzielić granice próbek siatki od fizycznych warunków brzegowych. Dla aktualnego stencilu wartość funkcji falowej jest zerowa na punktach pomocniczych o jeden krok poza skrajnymi próbkami.
- [ ] Dodać konwersje Å ↔ bohr i eV ↔ hartree tylko na granicy API; wewnątrz obliczeń używać jednostek atomowych.
- [ ] Ustalić strategię dla jądra trafiającego dokładnie w węzeł siatki: nie może powstać `Inf`/`NaN` w potencjale Coulomba. Na etapie prototypu można umieścić jądra między węzłami i mierzyć zależność wyniku od przesunięcia siatki; docelowo potrzebna jest metoda odporna na tę osobliwość.

Dla jednorodnej siatki prostokątnej:

$$
\Delta V = |h_x h_y h_z|,\qquad
\sum_{\mathbf r_g} |\psi_i(\mathbf r_g)|^2\,\Delta V \approx 1.
$$

`eigen_solver` zwraca wektor zwykle znormalizowany euklidesowo, `Σ_g |v_{ig}|² = 1`. Wartości fizyczne na węzłach to wtedy `ψ_i(r_g) = v_{ig}/√ΔV`. Ten przelicznik ma być wykonany w rdzeniu, zanim powstanie gęstość.

**Gotowe, gdy:** testy sprawdzają kształt danych, pozycje punktów, jednostki, normalizację fizyczną i poprawne zachowanie dla wadliwej geometrii.

### 2. Dokończyć niezależny solver jednoelektronowy

- [ ] Dodać ustawianą tolerancję, limit iteracji, kontrolę zbieżności i jawny raport residuum do `eigen_solver`. Obecna tolerancja `1e-4` jest zaszyta w wywołaniu KrylovKit.
- [ ] Dla każdego stanu zapisywać wartość własną, wektor, normę residuum i informację o zbieżności. Ustalić deterministyczne inicjowanie do testów.
- [ ] Sprawdzać ortogonalność stanów i normalizację z `ΔV`. Przy degeneracji porównywać podprzestrzeń/gęstość, ponieważ znaki i orientacja pojedynczych wektorów własnych mogą się zmieniać.
- [ ] Wprowadzić obiekt wyniku `OrbitalSet` z geometrią, poziomami energii, orbitalami, obsadzeniami i residuami. Liczba potrzebnych stanów niech zależy od liczby elektronów, a nie od stałej wartości domyślnej.
- [ ] Dodać obliczenie układów: pudełko, oscylator 3D, atom H, a potem H₂⁺ z dwoma jądrami.
- [ ] Badać osobno zbieżność względem `h` i wielkości pudełka; dla układów z jądrami również przesunięcie siatki względem jąder.

Jednoelektronowy hamiltonian molekularny i jego równanie własne:

$$
\hat H_1 = -\frac12\nabla^2 - \sum_A\frac{Z_A}{|\mathbf r-\mathbf R_A|},
\qquad \hat H_1\psi_i=\varepsilon_i\psi_i.
$$

Dla H₂⁺ energia całkowita przy ustalonej odległości jąder $R$ to energia najniższego elektronu **plus odpychanie jąder**: $E_{\mathrm{H_2^+}}(R)=\varepsilon_0(R)+1/R$. Dla oscylatora $V(r)=\frac12\omega^2 r^2$ poziom podstawowy w nieskończonej przestrzeni ma energię $E_0=3\omega/2$.

Norma residuum, którą raportujemy, może mieć postać:

$$
r_i=\frac{\|H v_i-\varepsilon_i v_i\|_2}
          {\max(1,\|H\|\,\|v_i\|_2)}.
$$

Normę operatora można zastąpić tańszą, jasno udokumentowaną skalą. Nie uznawać samego zakończenia iteracji za dowód poprawności stanu.

**Gotowe, gdy:** H i H₂⁺ dają stabilne stany podstawowe, a testy energii i residuów przechodzą dla kilku siatek.

### 3. Pierwsza wizualizacja prawdziwych orbitali molekularnych

Ten etap korzysta z `OrbitalSet` z etapu 2 i pozostaje w rozszerzeniu graficznym.

- [ ] Przeliczać wierzchołki `Meshing.isosurface` z indeksów tablicy na bohr: `r=(x_1+(i-1)h_x,\ y_1+(j-1)h_y,\ z_1+(k-1)h_z)`. Sprawdzić na polu o znanej izopowierzchni, czy biblioteka numeruje wierzchołki od 0 czy od 1.
- [ ] Pokazywać dodatni i ujemny płat rzeczywistego orbitalu różnymi kolorami; rysować też `|ψ|²` oraz pozycje jąder.
- [ ] Opisywać osie i odległości w bohr albo jawnie przeliczonych Å. Dopilnować poprawnych proporcji dla `hx ≠ hy ≠ hz`.
- [ ] Umożliwić wybór poziomu izopowierzchni oraz opcjonalnie poziomu zawierającego zadaną część prawdopodobieństwa, wyznaczaną z `Σ_g |ψ_g|² ΔV`.
- [ ] Zapewnić eksport obrazu i danych pola do standardowego formatu wolumetrycznego. Eksport danych ma być dostępny z rdzenia także bez wyświetlacza.
- [ ] Dodać testy matematycznego odwzorowania współrzędnych; test tworzenia figury może działać osobno od testów rdzenia.

**Gotowe, gdy:** orbitale H₂⁺ są umieszczone w tych samych współrzędnych co jądra i wygląd nie zmienia się błędnie po zmianie kroku jednej osi.

### 4. Solver Poissona dla potencjału Hartree

- [ ] Przyjmować gęstość `ρ_g` na `Grid3D` i zwracać potencjał `v_H` w tych samych węzłach.
- [ ] Na początek zaimplementować prosty solver z warunkiem brzegowym i dużym pudełkiem, ale mierzyć błąd po zmianie pudełka. Dla izolowanej cząsteczki docelowo dodać warunek asymptotyczny z ładunku/momentów multipolowych lub większą siatkę pomocniczą.
- [ ] Używać tej samej konwencji operatora Laplace’a i położenia brzegów co w części jednoelektronowej; dokumentować różnicę między warunkiem dla `ψ` a warunkiem dla `v_H`.
- [ ] Testować potencjał znanej, gładkiej gęstości oraz energię Hartree. Sprawdzać symetrię i zbieżność względem `h` oraz pudełka.
- [ ] Dopiero po poprawnym solverze wybrać przyspieszenie, np. metodę wielosiatkową lub solver iteracyjny ze skutecznym preconditionerem.

$$
-\nabla^2 v_H(\mathbf r)=4\pi\rho(\mathbf r),\qquad
E_H[\rho]=\frac12\int\rho(\mathbf r)v_H(\mathbf r)\,d^3r
       \approx\frac12\sum_g\rho_g(v_H)_g\Delta V.
$$

Przykład testowy: dla $\rho(r)=(\alpha/\pi)^{3/2}e^{-\alpha r^2}$ i warunku izolowanego układu $v_H(r)=\text{erf}(\sqrt{\alpha}\,r)/r$ dla $r>0$, a $v_H(0)=2\sqrt{\alpha/\pi}$. Błąd na brzegu pudełka należy analizować oddzielnie od błędu dyskretyzacji.

**Gotowe, gdy:** potencjał i energia testowej gęstości zbiegają do znanego wyniku, a błąd warunku brzegowego jest mierzony.

### 5. Składniki energii i prosty funkcjonał XC

- [ ] Dodać interfejs funkcjonału zwracający **razem** $E_{\rm xc}[\rho]$ i zgodny z nim $v_{\rm xc}(r)=\delta E_{\rm xc}/\delta\rho(r)$; osobne implementacje łatwo rozjeżdżają się numerycznie.
- [ ] Najpierw zaimplementować wymianę lokalną, potem dodać korelację LDA. Referencyjną implementację korelacji można porównać z Libxc.
- [ ] Obsłużyć `ρ=0` bez `NaN` i bez przypadkowej utraty znaku/typu. Przetestować pochodną energii numerycznie: zmiana $E_{\rm xc}$ po małym zaburzeniu gęstości powinna zgadzać się z $v_{\rm xc}$.
- [ ] Przechowywać osobno $T_s$, $E_{\rm ext}$, $E_H$, $E_{\rm xc}$ i $E_{nn}$ oraz ich sumę. Nie traktować sumy wartości własnych jako energii całkowitej.

Dla spinowo niespolaryzowanej wymiany LDA:

$$
E_x[\rho]=-C_x\int\rho^{4/3}(\mathbf r)\,d^3r,\qquad
C_x=\frac34\left(\frac3\pi\right)^{1/3},
$$

$$
v_x(\mathbf r)=\frac{\delta E_x}{\delta\rho(\mathbf r)}
=-\left(\frac3\pi\right)^{1/3}\rho(\mathbf r)^{1/3}.
$$

Wynik z samą wymianą opisywać jawnie jako **KS-LDA(x)**; wynik z konkretną parametryzacją korelacji jako **KS-LDA(XC)**. LDA dla jednego elektronu zawiera błąd samooddziaływania, więc dokładna energia H/H₂⁺ nie jest właściwym testem równości energii KS-LDA.

**Gotowe, gdy:** energia, potencjał XC i ich pochodna są zgodne w testach; raport energii pokazuje wszystkie składniki.

### 6. Pętla samouzgodnienia KS dla powłok zamkniętych

- [ ] Dodać obsadzenia `f_i∈{0,2}` dla par elektronów oraz inicjalną gęstość o całce `N_e`.
- [ ] Jeden krok SCF ma: z $\rho_k$ wyliczyć $v_H$ i $v_{\rm xc}$, zbudować efektywny hamiltonian, obliczyć potrzebne stany, znormalizować je, utworzyć $\rho_{\rm out}$ i energię.
- [ ] Zacząć od tłumionego mieszania `ρ_{k+1}=(1-α)ρ_k+αρ_{\rm out}`. Potem dodać Andersona/DIIS i ograniczenia zapobiegające ujemnej gęstości lub utracie liczby elektronów.
- [ ] Wymagać równocześnie małego błędu gęstości, stabilnej energii i małych residuów stanów własnych. Zwracać `converged=false` z diagnostyką po osiągnięciu `maxiter`, zamiast zwracać pozornie poprawną liczbę.
- [ ] Zwracać `SCFResult`: geometria, siatka, orbitale, obsadzenia, `ρ`, potencjały, energie, liczba iteracji, historia residuów i kryteria stopu.
- [ ] Porównać He i H₂ z referencją używającą **tej samej** definicji LDA, geometrii, jednostek i konwencji energii.

Równania do implementacji:

$$
\left[-\frac12\nabla^2+v_{\rm ext}(\mathbf r)
      +v_H[\rho](\mathbf r)+v_{\rm xc}[\rho](\mathbf r)\right]\psi_i
      =\varepsilon_i\psi_i,
$$

$$
\rho(\mathbf r)=\sum_i f_i|\psi_i(\mathbf r)|^2,\qquad
\sum_g\rho_g\Delta V \approx N_e,\qquad
v_{\rm ext}(\mathbf r)=-\sum_A\frac{Z_A}{|\mathbf r-\mathbf R_A|}.
$$

Dla układu bez smearing i z powłoką zamkniętą można obliczać energię całkowitą także ze wzoru kontrolnego:

$$ E_{\rm tot} = \sum_i f_i\varepsilon_i - E_H[\rho] + E_{\rm xc}[\rho] - \int \rho(\mathbf r)\,v_{\rm xc}(\mathbf r)\,d^3r + E_{nn}, 
$$

$$
E_{nn} = \sum_{A\lt B} \frac{Z_A Z_B} {\lvert \mathbf R_A-\mathbf R_B\rvert}
$$

Porównać tę postać z sumą składników $T_s+E_{\rm ext}+E_H+E_{\rm xc}+E_{nn}$, aby wykryć błąd podwójnego liczenia.

**Gotowe, gdy:** He i H₂ osiągają zadane kryteria, `∫ρ=N_e`, a obie drogi liczenia `E_{\rm tot}` zgadzają się w granicach błędu numerycznego.

### 7. Walidacja fizyczna i pierwszy wynik projektu

- [ ] Obliczyć krzywą `E_{\rm tot}(R)` dla H₂: kilka odległości jąder, ta sama metoda XC i kryteria zbieżności w każdym punkcie.
- [ ] Wyznaczyć minimum `R_{\rm eq}` oraz energię wiązania względem dwóch odseparowanych atomów H. Dla energii atomu H uwzględnić, że układ ma niesparowany elektron; jeśli solver nie obsługuje jeszcze spinu, w pierwszej wersji skupić się na krzywej i minimum H₂, a energię wiązania dodać po etapie spinowym.
- [ ] Zrobić tabelę zbieżności względem `h` i rozmiaru pudełka oddzielnie. Dla potencjału Coulomba dodać test przesunięcia siatki względem jąder.
- [ ] Porównać energię i gęstość z PySCF dla takiej samej geometrii i funkcjonału, po sprawdzeniu zbieżności bazy PySCF; nie porównywać orbitali punkt po punkcie bez uwzględnienia dowolnego znaku i obrotu w podprzestrzeni zdegenerowanej.
- [ ] Przygotować uruchamialny przykład bez grafiki zapisujący wynik i tabelę energii; oddzielny przykład graficzny odczytuje ten wynik.

**Gotowe, gdy:** użytkownik może odtworzyć krzywą H₂ jednym skryptem, a raport pokazuje błąd siatki, pudełka i SCF.

### 8. Skalowanie i późniejsze funkcje

Po uzyskaniu poprawnych wyników He/H₂:

- [ ] Wdrożyć macierzowo wolne działanie operatora kinetycznego na wektor; porównać czas i pamięć z obecną macierzą rzadką.
- [ ] Dodać wyższy rząd stencilu i testy zbieżności; nie zakładać automatycznie monotonicznej zbieżności energii metody różnic skończonych.
- [ ] Rozważyć preconditioner solvera wartości własnych i wydajniejszy solver Poissona, gdy profilowanie pokaże wąskie gardło.
- [ ] Dodać spinowo spolaryzowany UKS dla atomów H, rodników i układów z nieparzystą liczbą elektronów, a potem korelację LDA i GGA/PBE.
- [ ] Dodać odczyt XYZ, zapis/wznowienie `SCFResult`, moment dipolowy $\boldsymbol\mu=\sum_A Z_A\mathbf R_A-\int\mathbf r\rho(\mathbf r)\,d^3r$, siły i optymalizację geometrii. Siły wymagają osobnych testów wpływu przesunięcia jąder względem siatki.
- [ ] Większe cząsteczki, hybrydy i pseudopotencjały rozważać dopiero po walidacji prostszych etapów.

## Minimalny zestaw testów przekrojowych

1. **Matematyka siatki:** kształt, kolejność punktów, kroki osi, granice i `ΔV`.
2. **Solver:** pudełko 3D ($E_0=\frac{\pi^2}{2}(L_x^{-2}+L_y^{-2}+L_z^{-2})$), oscylator, residua, ortogonalność, normalizacja.
3. **Coulomb:** brak `Inf`/`NaN`, zbieżność H/H₂⁺ przy zmianie `h`, pudełka i przesunięcia siatki.
4. **Poisson:** rozkład Gaussa, energia Hartree, wpływ warunku brzegowego.
5. **KS:** zachowanie liczby elektronów, zgodność dwóch wzorów na energię, zgodność `E_{\rm xc}` z `v_{\rm xc}`, udokumentowana zbieżność SCF.
6. **Rozdział pakietu:** testy rdzenia bez zależności graficznych, osobny test renderowania oraz zgodność danych przekazywanych do rozszerzenia.

## Źródła techniczne do implementacji

- [DFTK: równania KS i SCF](https://docs.dftk.org/dev/guide/self_consistent_field/) — struktura pętli i mieszanie gęstości. DFTK koncentruje się na układach okresowych; tutaj warunki brzegowe trzeba dobrać dla izolowanych cząsteczek.
- [PySCF: molecular DFT](https://pyscf.org/user/dft.html) — definicje funkcjonałów i obliczenia referencyjne.
- [GPAW: solver Poissona](https://gpaw.readthedocs.io/documentation/poisson.html) — wpływ warunków brzegowych i poprawki multipolowe.
- [GPAW: reprezentacja na siatce](https://gpaw.readthedocs.io/tutorialsexercises/basics/representation/representation.html) — zbieżność metod różnic skończonych.
- [Libxc: lista funkcjonałów](https://libxc.gitlab.io/functionals/) — wybór i porównanie parametryzacji XC.
- [Julia Pkg: słabe zależności i rozszerzenia](https://pkgdocs.julialang.org/v1/creating-packages/#Weak-dependencies) — sposób technicznego oddzielenia wykresów od rdzenia.
