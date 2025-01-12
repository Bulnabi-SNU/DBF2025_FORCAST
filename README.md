# DBF2025_SizingCode

Analysis team code used for the AIAA DBF 2025 competition.

### VSP_analysis 사용법

* [초기설정 가이드](https://www.notion.so/Setup-Guide-e91707f7552e4c849eb1648c549a9c9c)를 따라하기

* vsp_analysis.py에서 airfoil data path를 적절히 수정

```bash
s9027_path = r"C:\Path\to\DBF2025_SizingCode\VSP_analysis\s9027.dat"
naca0008_path = r"C:\Path\to\DBF2025_SizingCode\VSP_analysis\naca0008.dat"
```

* Anaylsis parameter 설정

```bash
vsp_file = os.path.join(vsp3_dir, "Mothership.vsp3")
alpha_start = -2.0   # Starting angle of attack (degrees)
alpha_end = 2.0      # Ending angle of attack (degrees)
alpha_step = 1.0     # Step size (degrees)
Re = 380000          # Reynolds number
Mach = 0             # Mach number (subsonic)
```