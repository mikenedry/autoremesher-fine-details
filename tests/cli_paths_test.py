"""Standalone Unicode path regressions for the actual CLI executable."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from cli_test_mesh import grid

ENGINE = Path(sys.argv.pop(1)).resolve()


class CLIPaths(unittest.TestCase):
    def invoke(self, source, output, report=None, cwd=None):
        args = [str(ENGINE), '--input', str(source), '--output', str(output),
                '--quads', '64', '--no-cloth-folds', '--quiet']
        if report is not None:
            args += ['--report', str(report)]
        return subprocess.run(args, cwd=cwd, capture_output=True, encoding='utf-8',
                              errors='replace', timeout=120)

    def test_empty_input_reaches_validation_in_unicode_directory(self):
        with tempfile.TemporaryDirectory(prefix='autoremesher paths ') as root:
            for name in ('ascii', 'Jos\u00e9', '\u6d4b\u8bd5', '\U0001f680 space'):
                with self.subTest(name=name):
                    folder = Path(root)/name
                    folder.mkdir()
                    source = folder/'scan.obj'
                    source.write_text('# Empty OBJ\n', encoding='ascii')
                    result = self.invoke(source, folder/'output.obj')
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertIn('Empty input mesh', result.stderr)
                    self.assertFalse((folder/'output.obj').exists())

    def test_unicode_input_output_and_report_roundtrip(self):
        with tempfile.TemporaryDirectory(prefix='autoremesher paths ') as root:
            for boundary in ('input', 'output', 'report', 'all'):
                with self.subTest(boundary=boundary):
                    folder = Path(root)/boundary
                    folder.mkdir()
                    paths = {}
                    for kind, suffix in (('input', '.obj'), ('output', '.obj'), ('report', '.json')):
                        name = '\u6d4b\u8bd5 Jos\u00e9 \U0001f680' if boundary in (kind, 'all') else 'ascii'
                        parent = folder/(kind+' '+name)
                        parent.mkdir()
                        paths[kind] = parent/(name+suffix)
                    source = grid()
                    paths['input'].write_text(source, encoding='ascii')
                    result = self.invoke(paths['input'], paths['output'], paths['report'], cwd=folder)
                    self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                    report = json.loads(paths['report'].read_text(encoding='utf-8'))
                    self.assertEqual(report, json.loads(result.stdout.strip().splitlines()[-1]))
                    obj = paths['output'].read_text(encoding='ascii').splitlines()
                    self.assertEqual(report['vertices'], sum(line.startswith('v ') for line in obj))
                    self.assertEqual(report['faces'], sum(line.startswith('f ') for line in obj))
                    self.assertGreater(report['quads'], 0)
                    self.assertEqual(paths['input'].read_text(encoding='ascii'), source)


if __name__ == '__main__':
    unittest.main(verbosity=2)
