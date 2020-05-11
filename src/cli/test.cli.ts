import CliHelper, { CliListener } from "interactive-cli-helper";
import { Martinizer } from "../Builders/Martinizer";
import path from "path";
import { promises as FsPromise } from 'fs';
import { LIPIDS_ROOT_DIR } from "../constants";
import { Lipid } from "../Entities/entities";
import { Database } from "../Entities/CouchHelper";
import RadiusDatabase from "../Entities/RadiusDatabase";
import MembraneBuilder from "../Builders/MembraneBuilder";
import TmpDirHelper from "../TmpDirHelper/TmpDirHelper";

const TEST_CLI = new CliListener(CliHelper.formatHelp('test', {
  map: 'Create a map with the web server with the specified file',
  ccmap: 'Create a map with the Python CCMap package with the specified file',
  'elastic-bonds': 'Generate the elastic bonds definitions with the given folder. It must contain a TOP file, PDB file and ITP(s)',
  'go-sites': 'Generate the elastic bonds definitions with the given folder. It must contain a ITPs files describing a Go model',
  'add-lipid': 'Add a lipid, followed by a lipid filename (search in {workdir}/lipids/2_2).',
  'lipid-prefixes': 'Get lipid directory prefix name by force field.',
  'auto-import-lipid <forceField>': 'Import all the itp files inside a lipid directory found by its related force field.',
  'dry-run': 'Run a Martinize + Go + INSANE + GROMACS test in order to test features.',
}));

TEST_CLI.addSubListener('map', rest => Martinizer.getMap(rest));

TEST_CLI.addSubListener('ccmap', rest => Martinizer.getCcMap(rest));

TEST_CLI.addSubListener('elastic-bonds', async rest => {
  const itps: string[] = [];
  let top_file: string = "", pdb_file: string = "";

  rest = path.resolve(rest) + "/";
  console.log("Given path:", rest);

  for (const element of await FsPromise.readdir(rest)) {
    const basename = path.basename(element);
    if (basename.endsWith('.pdb')) {
      pdb_file = rest + basename;
    }
    else if (basename.endsWith('.top')) {
      top_file = rest + basename;
    }
    else if (basename.endsWith('.itp')) {
      itps.push(rest + basename);
    }
  }

  if (!top_file || !pdb_file || !itps.length) {
    return "Given folder must have a ITP file, TOP file and PDB file.";
  }

  const relations = await Martinizer.computeElasticNetworkBounds(top_file, itps);

  await FsPromise.writeFile(rest + 'relations.json', JSON.stringify(relations, null, 2));

  return "Relations has been written to '" + (rest + 'relations.json') + "'.";
});

TEST_CLI.addSubListener('go-sites', async rest => {
  const itps: string[] = [];
  let top_file: string = "";

  rest = path.resolve(rest) + "/";
  console.log("Given path:", rest);

  for (const element of await FsPromise.readdir(rest)) {
    const basename = path.basename(element);
    if (basename.endsWith('.top')) {
      top_file = rest + basename;
    }
    else if (basename.endsWith('.itp')) {
      itps.push(rest + basename);
    }
  }

  if (!itps.length) {
    return "Given folder must have ITP files.";
  }

  const relations = await Martinizer.__UNSAFEcomputeGoModelBounds(top_file, itps);

  await FsPromise.writeFile(rest + 'relations.json', JSON.stringify(relations, null, 2));

  return "Relations has been written to '" + (rest + 'relations.json') + "'.";
});

TEST_CLI.addSubListener('add-lipid', async rest => {
  const name = LIPIDS_ROOT_DIR + "2_2/" + rest + '.itp';

  const lipid: Lipid = {
    id: '',
    name: rest,
    itp: await FsPromise.readFile(name, 'utf-8'),
  };

  await Database.lipid.add(lipid, 'martini22');

  return lipid;
});

TEST_CLI.addSubListener('lipid-prefixes', () => {
  return CliHelper.formatHelp("    force field  prefix", RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION); 
});

TEST_CLI.addSubListener('auto-import-lipid', async force_field => {
  const prefix = RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[force_field];

  if (!prefix) {
    return "Unknown force field.";
  }

  const dir = LIPIDS_ROOT_DIR + prefix + "/";

  for (const itp of await FsPromise.readdir(dir)) {
    if (!itp.endsWith('.itp')) {
      continue;
    }

    const lipid: Lipid = {
      id: '',
      name: itp.split('.itp')[0].toLocaleUpperCase(),
      itp: await FsPromise.readFile(dir + itp, 'utf-8'),
    };

    await Database.lipid.add(lipid, force_field);
  }
});

export default TEST_CLI;

const KWALP_TEST_PDB = `
ATOM      1  N   GLY A   1       5.974 -11.703   1.119  1.00101.82           N
ATOM      2  CA  GLY A   1       4.991 -12.260   0.070  1.00101.82           C
ATOM      3  C   GLY A   1       3.687 -11.417  -0.173  1.00101.82           C
ATOM      4  O   GLY A   1       3.770 -10.226  -0.472  1.00101.82           O
ATOM      8  N   LYS A   2       2.508 -12.065  -0.050  1.00 95.66           N
ATOM      9  CA  LYS A   2       1.269 -11.377  -0.278  1.00 95.66           C
ATOM     10  CB  LYS A   2      -0.016 -12.201  -0.074  1.00 95.66           C
ATOM     11  CG  LYS A   2      -0.367 -13.207  -1.164  1.00 95.66           C
ATOM     12  CD  LYS A   2      -1.640 -13.987  -0.831  1.00 95.66           C
ATOM     13  CE  LYS A   2      -2.641 -14.025  -1.987  1.00 95.66           C
ATOM     14  NZ  LYS A   2      -2.758 -12.682  -2.602  1.00 95.66           N
ATOM     15  C   LYS A   2       1.128 -10.301   0.736  1.00 95.66           C
ATOM     16  O   LYS A   2       0.661  -9.212   0.413  1.00 95.66           O
ATOM     17  N   LYS A   3       1.513 -10.604   1.987  1.00 97.39           N
ATOM     18  CA  LYS A   3       1.361  -9.717   3.107  1.00 97.39           C
ATOM     19  CB  LYS A   3       2.023 -10.320   4.359  1.00 97.39           C
ATOM     20  CG  LYS A   3       1.376 -11.585   4.927  1.00 97.39           C
ATOM     21  CD  LYS A   3       0.106 -11.355   5.747  1.00 97.39           C
ATOM     22  CE  LYS A   3       0.404 -10.898   7.178  1.00 97.39           C
ATOM     23  NZ  LYS A   3      -0.602 -11.456   8.106  1.00 97.39           N
ATOM     24  C   LYS A   3       2.145  -8.490   2.851  1.00 97.39           C
ATOM     25  O   LYS A   3       1.660  -7.382   3.068  1.00 97.39           O
ATOM     26  N   LYS A   4       3.391  -8.678   2.381  1.00103.17           N
ATOM     27  CA  LYS A   4       4.289  -7.585   2.185  1.00103.17           C
ATOM     28  CB  LYS A   4       5.702  -8.026   1.771  1.00103.17           C
ATOM     29  CG  LYS A   4       5.825  -8.393   0.294  1.00103.17           C
ATOM     30  CD  LYS A   4       7.272  -8.535  -0.185  1.00103.17           C
ATOM     31  CE  LYS A   4       8.219  -7.471   0.371  1.00103.17           C
ATOM     32  NZ  LYS A   4       8.795  -7.939   1.651  1.00103.17           N
ATOM     33  C   LYS A   4       3.765  -6.686   1.118  1.00103.17           C
ATOM     34  O   LYS A   4       3.845  -5.464   1.250  1.00103.17           O
ATOM     35  N   LEU A   5       3.221  -7.269   0.027  1.00129.32           N
ATOM     36  CA  LEU A   5       2.767  -6.450  -1.060  1.00129.32           C
ATOM     37  CB  LEU A   5       2.277  -7.243  -2.289  1.00129.32           C
ATOM     38  CG  LEU A   5       1.626  -6.349  -3.369  1.00129.32           C
ATOM     39  CD1 LEU A   5       2.574  -5.229  -3.821  1.00129.32           C
ATOM     40  CD2 LEU A   5       1.104  -7.183  -4.552  1.00129.32           C
ATOM     41  C   LEU A   5       1.641  -5.587  -0.624  1.00129.32           C
ATOM     42  O   LEU A   5       1.587  -4.407  -0.970  1.00129.32           O
ATOM     43  N   ALA A   6       0.715  -6.172   0.155  1.00 31.40           N
ATOM     44  CA  ALA A   6      -0.450  -5.465   0.596  1.00 31.40           C
ATOM     45  CB  ALA A   6      -1.399  -6.374   1.391  1.00 31.40           C
ATOM     46  C   ALA A   6      -0.056  -4.313   1.460  1.00 31.40           C
ATOM     47  O   ALA A   6      -0.631  -3.225   1.350  1.00 31.40           O
ATOM     48  N   LEU A   7       0.928  -4.524   2.356  1.00129.53           N
ATOM     49  CA  LEU A   7       1.310  -3.471   3.258  1.00129.53           C
ATOM     50  CB  LEU A   7       2.282  -3.911   4.359  1.00129.53           C
ATOM     51  CG  LEU A   7       2.525  -2.782   5.374  1.00129.53           C
ATOM     52  CD1 LEU A   7       1.194  -2.328   5.999  1.00129.53           C
ATOM     53  CD2 LEU A   7       3.564  -3.184   6.433  1.00129.53           C
ATOM     54  C   LEU A   7       1.916  -2.328   2.499  1.00129.53           C
ATOM     55  O   LEU A   7       1.683  -1.168   2.841  1.00129.53           O
ATOM     56  N   ALA A   8       2.731  -2.627   1.464  1.00 30.38           N
ATOM     57  CA  ALA A   8       3.367  -1.579   0.709  1.00 30.38           C
ATOM     58  CB  ALA A   8       4.308  -2.108  -0.390  1.00 30.38           C
ATOM     59  C   ALA A   8       2.335  -0.742   0.029  1.00 30.38           C
ATOM     60  O   ALA A   8       2.452   0.482  -0.008  1.00 30.38           O
ATOM     61  N   LEU A   9       1.291  -1.393  -0.525  1.00127.21           N
ATOM     62  CA  LEU A   9       0.267  -0.690  -1.254  1.00127.21           C
ATOM     63  CB  LEU A   9      -0.762  -1.685  -1.846  1.00127.21           C
ATOM     64  CG  LEU A   9      -2.021  -1.087  -2.512  1.00127.21           C
ATOM     65  CD1 LEU A   9      -3.036  -0.563  -1.478  1.00127.21           C
ATOM     66  CD2 LEU A   9      -1.662  -0.033  -3.572  1.00127.21           C
ATOM     67  C   LEU A   9      -0.446   0.247  -0.340  1.00127.21           C
ATOM     68  O   LEU A   9      -0.736   1.381  -0.723  1.00127.21           O
ATOM     69  N   ALA A  10      -0.753  -0.213   0.887  1.00 34.55           N
ATOM     70  CA  ALA A  10      -1.485   0.589   1.828  1.00 34.55           C
ATOM     71  CB  ALA A  10      -1.790  -0.163   3.130  1.00 34.55           C
ATOM     72  C   ALA A  10      -0.694   1.807   2.192  1.00 34.55           C
ATOM     73  O   ALA A  10      -1.263   2.892   2.318  1.00 34.55           O
ATOM     74  N   LEU A  11       0.633   1.653   2.397  1.00 48.49           N
ATOM     75  CA  LEU A  11       1.429   2.779   2.795  1.00 48.49           C
ATOM     76  CB  LEU A  11       2.881   2.442   3.181  1.00 48.49           C
ATOM     77  CG  LEU A  11       2.984   1.629   4.487  1.00 48.49           C
ATOM     78  CD1 LEU A  11       4.441   1.510   4.968  1.00 48.49           C
ATOM     79  CD2 LEU A  11       2.043   2.187   5.567  1.00 48.49           C
ATOM     80  C   LEU A  11       1.456   3.812   1.713  1.00 48.49           C
ATOM     81  O   LEU A  11       1.388   5.011   1.990  1.00 48.49           O
ATOM     82  N   ALA A  12       1.566   3.368   0.446  1.00 31.04           N
ATOM     83  CA  ALA A  12       1.633   4.308  -0.639  1.00 31.04           C
ATOM     84  CB  ALA A  12       1.809   3.627  -2.005  1.00 31.04           C
ATOM     85  C   ALA A  12       0.370   5.099  -0.716  1.00 31.04           C
ATOM     86  O   ALA A  12       0.411   6.310  -0.943  1.00 31.04           O
ATOM     87  N   LEU A  13      -0.787   4.427  -0.536  1.00 47.31           N
ATOM     88  CA  LEU A  13      -2.058   5.096  -0.640  1.00 47.31           C
ATOM     89  CB  LEU A  13      -3.263   4.147  -0.534  1.00 47.31           C
ATOM     90  CG  LEU A  13      -3.464   3.316  -1.809  1.00 47.31           C
ATOM     91  CD1 LEU A  13      -4.726   2.444  -1.730  1.00 47.31           C
ATOM     92  CD2 LEU A  13      -3.454   4.227  -3.050  1.00 47.31           C
ATOM     93  C   LEU A  13      -2.186   6.132   0.427  1.00 47.31           C
ATOM     94  O   LEU A  13      -2.716   7.217   0.171  1.00 47.31           O
ATOM     95  N   ALA A  14      -1.733   5.822   1.660  1.00 29.40           N
ATOM     96  CA  ALA A  14      -1.855   6.760   2.745  1.00 29.40           C
ATOM     97  CB  ALA A  14      -1.362   6.187   4.087  1.00 29.40           C
ATOM     98  C   ALA A  14      -1.053   7.987   2.456  1.00 29.40           C
ATOM     99  O   ALA A  14      -1.503   9.102   2.729  1.00 29.40           O
ATOM    100  N   LEU A  15       0.165   7.811   1.900  1.00 42.55           N
ATOM    101  CA  LEU A  15       1.006   8.947   1.652  1.00 42.55           C
ATOM    102  CB  LEU A  15       2.404   8.577   1.106  1.00 42.55           C
ATOM    103  CG  LEU A  15       3.267   7.716   2.054  1.00 42.55           C
ATOM    104  CD1 LEU A  15       4.719   7.628   1.551  1.00 42.55           C
ATOM    105  CD2 LEU A  15       3.169   8.188   3.513  1.00 42.55           C
ATOM    106  C   LEU A  15       0.366   9.845   0.642  1.00 42.55           C
ATOM    107  O   LEU A  15       0.391  11.067   0.789  1.00 42.55           O
ATOM    108  N   ALA A  16      -0.224   9.250  -0.413  1.00 28.83           N
ATOM    109  CA  ALA A  16      -0.816  10.036  -1.460  1.00 28.83           C
ATOM    110  CB  ALA A  16      -1.386   9.173  -2.596  1.00 28.83           C
ATOM    111  C   ALA A  16      -1.953  10.843  -0.918  1.00 28.83           C
ATOM    112  O   ALA A  16      -2.126  12.004  -1.294  1.00 28.83           O
ATOM    113  N   LEU A  17      -2.776  10.243  -0.034  1.00 40.59           N
ATOM    114  CA  LEU A  17      -3.929  10.939   0.478  1.00 40.59           C
ATOM    115  CB  LEU A  17      -4.815  10.044   1.364  1.00 40.59           C
ATOM    116  CG  LEU A  17      -5.457   8.881   0.586  1.00 40.59           C
ATOM    117  CD1 LEU A  17      -6.466   8.114   1.455  1.00 40.59           C
ATOM    118  CD2 LEU A  17      -6.062   9.369  -0.741  1.00 40.59           C
ATOM    119  C   LEU A  17      -3.498  12.116   1.288  1.00 40.59           C
ATOM    120  O   LEU A  17      -4.112  13.184   1.225  1.00 40.59           O
ATOM    121  N   TRP A  18      -2.430  11.941   2.088  1.00 96.93           N
ATOM    122  CA  TRP A  18      -2.023  13.023   2.927  1.00 96.93           C
ATOM    123  CB  TRP A  18      -0.942  12.673   3.961  1.00 96.93           C
ATOM    124  CG  TRP A  18      -0.750  13.800   4.946  1.00 96.93           C
ATOM    125  CD2 TRP A  18       0.192  14.863   4.755  1.00 96.93           C
ATOM    126  CD1 TRP A  18      -1.392  14.058   6.124  1.00 96.93           C
ATOM    127  NE1 TRP A  18      -0.908  15.220   6.678  1.00 96.93           N
ATOM    128  CE2 TRP A  18       0.068  15.726   5.844  1.00 96.93           C
ATOM    129  CE3 TRP A  18       1.082  15.099   3.749  1.00 96.93           C
ATOM    130  CZ2 TRP A  18       0.844  16.847   5.940  1.00 96.93           C
ATOM    131  CZ3 TRP A  18       1.866  16.226   3.850  1.00 96.93           C
ATOM    132  CH2 TRP A  18       1.745  17.082   4.925  1.00 96.93           C
ATOM    133  C   TRP A  18      -1.545  14.172   2.093  1.00 96.93           C
ATOM    134  O   TRP A  18      -1.798  15.331   2.418  1.00 96.93           O
ATOM    135  N   TRP A  19      -0.810  13.876   1.006  1.00 57.06           N
ATOM    136  CA  TRP A  19      -0.285  14.924   0.173  1.00 57.06           C
ATOM    137  CB  TRP A  19       0.611  14.382  -0.947  1.00 57.06           C
ATOM    138  CG  TRP A  19       1.893  13.813  -0.403  1.00 57.06           C
ATOM    139  CD2 TRP A  19       3.061  13.557  -1.194  1.00 57.06           C
ATOM    140  CD1 TRP A  19       2.210  13.468   0.878  1.00 57.06           C
ATOM    141  NE1 TRP A  19       3.504  13.013   0.934  1.00 57.06           N
ATOM    142  CE2 TRP A  19       4.041  13.063  -0.334  1.00 57.06           C
ATOM    143  CE3 TRP A  19       3.291  13.724  -2.527  1.00 57.06           C
ATOM    144  CZ2 TRP A  19       5.275  12.723  -0.801  1.00 57.06           C
ATOM    145  CZ3 TRP A  19       4.542  13.380  -2.992  1.00 57.06           C
ATOM    146  CH2 TRP A  19       5.514  12.893  -2.144  1.00 57.06           C
ATOM    147  C   TRP A  19      -1.410  15.700  -0.457  1.00 57.06           C
ATOM    148  O   TRP A  19      -1.367  16.930  -0.509  1.00 57.06           O
ATOM    149  N   ALA A  20      -2.442  14.993  -0.956  1.00 88.42           N
ATOM    150  CA  ALA A  20      -3.533  15.653  -1.624  1.00 88.42           C
ATOM    151  CB  ALA A  20      -4.601  14.677  -2.144  1.00 88.42           C
ATOM    152  C   ALA A  20      -4.228  16.614  -0.664  1.00 88.42           C
ATOM    153  O   ALA A  20      -4.392  17.783  -1.052  1.00 88.42           O
TER     156      ALA A  20
END
`;

TEST_CLI.addSubListener('dry-run', async () => {
  // Martinize the KWALP (with Go sites, use GROMACS)
  const pdb_dir = await TmpDirHelper.get();

  // Save the pdb
  await FsPromise.writeFile(pdb_dir + '/test.pdb', KWALP_TEST_PDB);

  const { pdb, itps, top } = await Martinizer.run({
    input: pdb_dir + '/test.pdb',
    ff: 'martini304',
    position: 'backbone',
    use_go_virtual_sites: true,
  });

  // Insert in INSANE
  const insane = await MembraneBuilder.run({
    force_field: 'martini304',
    molecule_pdb: pdb,
    molecule_itps: itps,
    molecule_top: top,
    lipids: [['DLPC', 4], ['DPPC', 1]],
    upper_leaflet: [['DLPC', 1]],
    settings: {
      box: [7, 7, 9],
    }
  });

  await TmpDirHelper.revoke(pdb_dir);

  const dir = path.dirname(insane.top);

  return `Run end successfully, you can inspect out files in ${dir}.`;
});
