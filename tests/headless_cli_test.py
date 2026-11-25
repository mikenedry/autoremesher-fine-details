"""Exercise the independently built OBJ driver, including real geometry IO."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ENGINE = Path(sys.argv.pop(1)).resolve()


class HeadlessCLI(unittest.TestCase):
    def invoke(self, *args, cwd=None, timeout=120):
        return subprocess.run([str(ENGINE), *map(str, args)], cwd=cwd,
                              capture_output=True, text=True, timeout=timeout)

    def test_version_identifies_local_singularity_and_legacy_protocol(self):
        result = self.invoke('--version')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r'^autoremesher-fine-details (?:[0-9a-f]{40}(?:\+dirty)?|unknown) rebuilt-cli-2 local-singularity-1\n$')

    def test_help_describes_cloth_fold_default(self):
        result = self.invoke('--help')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--cloth-folds (default)', result.stdout)
        self.assertIn('--no-cloth-folds', result.stdout)

    def test_invalid_options_fail_before_creating_outputs(self):
        with tempfile.TemporaryDirectory(prefix='autoremesher cli validation ') as directory:
            output = Path(directory) / 'output.obj'
            for option, value in [('--quads', '0'), ('--quads', '-1'),
                                  ('--quads', '1.5'), ('--quads', 'nan'),
                                  ('--quads', 'inf'), ('--quads', '100000001'),
                                  ('--adaptivity', '-1'), ('--adaptivity', '2'),
                                  ('--model-type', 'unknown')]:
                with self.subTest(option=option, value=value):
                    result = self.invoke('--input', 'missing.obj', '--output', output,
                                         option, value, cwd=directory)
                    self.assertEqual(result.returncode, 2)
                    self.assertIn('AutoRemesherEngine:', result.stderr)
                    self.assertNotIn('OBJ load failed', result.stderr)
                    self.assertFalse(output.exists())

    def test_same_input_output_is_refused_without_modifying_source(self):
        with tempfile.TemporaryDirectory(prefix='autoremesher cli paths ') as directory:
            source = Path(directory) / 'source.obj'
            source.write_text('v 0 0 0\n')
            before = source.read_bytes()
            result = self.invoke('--input', source, '--output', source)
            self.assertEqual(result.returncode, 2)
            self.assertIn('Distinct --input and --output', result.stderr)
            self.assertEqual(source.read_bytes(), before)

    def test_open_grid_roundtrip_reports_actual_output(self):
        with tempfile.TemporaryDirectory(prefix='autoremesher OBJ é ') as directory:
            source = Path(directory) / 'input grid.obj'
            output = Path(directory) / 'output grid.obj'
            report = Path(directory) / 'report.json'
            side = 12
            lines = []
            for y in range(side + 1):
                for x in range(side + 1):
                    lines.append(f'v {x/side:.17g} {y/side:.17g} {0.06*math.sin(2*math.pi*x/side):.17g}')
            for y in range(side):
                for x in range(side):
                    a = y*(side+1)+x+1
                    lines += [f'f {a} {a+1} {a+side+2}', f'f {a} {a+side+2} {a+side+1}']
            source.write_text('\n'.join(lines)+'\n')
            original = hashlib.sha256(source.read_bytes()).hexdigest()
            result = self.invoke('--input', source, '--output', output, '--report', report,
                                 '--quads', '64', '--no-cloth-folds', '--quiet', cwd=directory)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            measured = json.loads(report.read_text())
            self.assertEqual(measured, json.loads(result.stdout.strip().splitlines()[-1]))
            self.assertFalse(measured['cloth_folds'])
            self.assertNotIn('PROGRESS ', result.stdout)
            vertices, faces = [], []
            for line in output.read_text().splitlines():
                fields = line.split()
                if fields[0] == 'v':
                    vertices.append([float(v) for v in fields[1:]])
                elif fields[0] == 'f':
                    faces.append([int(v)-1 for v in fields[1:]])
            self.assertTrue(vertices)
            self.assertTrue(faces)
            self.assertTrue(all(len(v)==3 and all(math.isfinite(x) for x in v) for v in vertices))
            self.assertTrue(all(len(f)>=3 and all(0<=i<len(vertices) for i in f) for f in faces))
            self.assertEqual(measured['vertices'], len(vertices))
            self.assertEqual(measured['faces'], len(faces))
            self.assertEqual(measured['quads'], sum(len(f)==4 for f in faces))
            self.assertGreater(measured['quads'], 0)
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), original)


if __name__ == '__main__':
    unittest.main(verbosity=2)
