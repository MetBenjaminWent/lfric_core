##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict, Optional, Tuple

from pytest import MonkeyPatch, fixture, mark, raises

from ..pkg_config import LinkType, Package, PackageException


@fixture
def system_path(tmp_path: Path) -> Path:
    system_path = tmp_path / 'usr' / 'local'
    system_path.mkdir(parents=True)
    return system_path


@fixture
def system_pkg_path(system_path: Path, monkeypatch: MonkeyPatch) -> Path:
    pkg_path = system_path / 'lib' / 'pkg_config'
    pkg_path.mkdir(parents=True)
    monkeypatch.setenv('PKG_CONFIG_LIBDIR', str(pkg_path))
    return pkg_path


@fixture
def user_path(tmp_path: Path) -> Path:
    user_path = tmp_path / 'opt' / 'special'
    user_path.mkdir(parents=True)
    return user_path


@fixture
def user_pkg_path(user_path: Path, monkeypatch: MonkeyPatch) -> Path:
    pkg_path = user_path / 'lib' / 'pkg_config'
    pkg_path.mkdir(parents=True)
    monkeypatch.setenv('PKG_CONFIG_PATH', str(pkg_path))
    return pkg_path


class TestPackage:
    @mark.parametrize('name, expected', [
        ("single", {
            'version': (1, 0, 1),
            'compile_args': ('-Iusr/local/include/simple',),
            'link_args': ('-Lusr/local/lib/simple', '-lsimple')
        }),
        ("pot_hole", {
            'version': (3, 2, 1),
            'compile_args': ('-Iusr/local/include/pothole',),
            'link_args': ('-Lusr/local/lib/pothole', '-lpothole')
        }),
        ("CamelCase", {
            'version': (1, 2, 3),
            'compile_args': ('-Iopt/special/include/camel',),
            'link_args': ('-Lopt/special/lib/camel', '-lcamel')
        })
    ])
    def test_constructor_name(self, name: str, expected: Dict[str, Any],
                              system_path: Path, system_pkg_path: Path,
                              user_path: Path, user_pkg_path: Path,
                              tmp_path: Path, monkeypatch: MonkeyPatch):
        """
        Checks some likely library names.

        Libraries on both system and user paths.
        """
        relative_path = system_path.relative_to(tmp_path)
        (system_pkg_path / 'single.pc').write_text(
            dedent(
                f"""
                Name: single
                Version: 1.0.1
                Description: Single word, lower case.
                URL: http://example.com/single
                Cflags: -I{relative_path}/include/simple
                Cflags.private: -I{relative_path}/include/static
                Libs: -L{relative_path}/lib/simple -lsimple
                Libs.private: -L{relative_path}/lib/static -l static
                """
            )
        )

        relative_path = system_path.relative_to(tmp_path)
        (system_pkg_path / 'pot_hole.pc').write_text(
            dedent(
                f"""
                Name: pot_hole
                Version: 3.2.1
                Description: Two words, unlerine separated.
                URL: http://example.com/pothole
                Cflags: -I{relative_path}/include/pothole
                Cflags.private: -I{relative_path}/include/static
                Libs: -L{relative_path}/lib/pothole -lpothole
                Libs.private: -L{relative_path}/lib/static -l static
                """
            )
        )

        relative_path = user_path.relative_to(tmp_path)
        (user_pkg_path / 'CamelCase.pc').write_text(
            dedent(
                f"""
                Name: CamelCase
                Version: 1.2.3
                Description: Two words, initial capital.
                URL: http://example.com/CamelCase
                Cflags: -I{relative_path}/include/camel
                Cflags.private: -I{relative_path}/include/static
                Libs: -L{relative_path}/lib/camel -lcamel
                Libs.private: -L{relative_path}/lib/static -l static
                """
            )
        )

        test_unit = Package(name)
        assert test_unit.name == name
        assert test_unit.version == expected['version']
        assert test_unit.compile_arguments == expected['compile_args']
        assert test_unit.link_arguments == expected['link_args']

    @mark.parametrize('version, expected', [
        ('', {
            'version': tuple()
        }),
        ('1', {
            'version': (1,)
        }),
        ('1.2', {
            'version': (1, 2)
        }),
        ('1.2.3', {
            'version': (1, 2, 3)
        }),
        ('1.2.dev1', {
            'version': (1, 2, 'dev1')
        }),
        ('2.3.dev.2', {
            'version': (2, 3, 'dev', 2)
        })
    ])
    def test_constructor_version(self, version: str,
                                 expected: Dict[str, Any],
                                 user_path: Path, user_pkg_path: Path,
                                 tmp_path: Path):
        """
        Checks various version number formats.

        Libraries on user paths.
        """
        (user_pkg_path / 'test.pc').write_text(
            dedent(
                f"""
                Name: test
                Version: {version}
                Description: Many versions.
                URL: http://example.com/version
                Cflags: -I{user_path.relative_to(tmp_path)}/include
                Cflags.private: -I{user_path.relative_to(tmp_path)}/include
                Libs: -L{user_path.relative_to(tmp_path)}/lib -lversion
                Libs.private: -L{user_path.relative_to(tmp_path)}/lib -l static
                """
            )
        )

        test_unit = Package('test')
        assert test_unit.name == 'test'
        assert test_unit.version == expected['version']
        assert test_unit.compile_arguments == ('-Iopt/special/include',)
        assert test_unit.link_arguments == ('-Lopt/special/lib', '-lversion')

    @mark.parametrize('arg_str, expected', [
        ('', tuple()),
        ('-I/usr/local/lib/special', ('-I/usr/local/lib/special',)),
        ('-I/usr/local/lib/special -I/usr/local/lib/other',
         ('-I/usr/local/lib/special', '-I/usr/local/lib/other')),
        ('-I /usr/local/special', ('-I/usr/local/special',))
    ])
    def test_constructor_compile_arguments(self, arg_str: str,
                                           expected: Tuple[str],
                                           user_path: Path,
                                           user_pkg_path: Path,
                                           tmp_path: Path):
        """
        Checks compiler path arguments are compressed.

        Libraries on user paths.
        """
        (user_pkg_path / 'test.pc').write_text(
            dedent(
                f"""
                Name: test
                Version: 1.0.0
                Description: Many versions.
                URL: http://example.com/version
                Cflags: {arg_str}
                Cflags.private:
                Libs: -L{user_path.relative_to(tmp_path)}/lib -lversion
                Libs.private: -L{user_path.relative_to(tmp_path)}/lib -l static
                """
            )
        )

        test_unit = Package('test')
        assert test_unit.name == 'test'
        assert test_unit.version == (1, 0, 0)
        assert test_unit.compile_arguments == expected
        assert test_unit.link_arguments == ('-Lopt/special/lib', '-lversion')

    @mark.parametrize('arg_str, expected', [
        ('', tuple()),
        ('-L/usr/local/lib/special -lspecial',
         ('-L/usr/local/lib/special', '-lspecial')),
        ('-L/usr/local/lib/special -lspecial -L/usr/local/lib/other -lother',
         ('-L/usr/local/lib/special', '-lspecial',
          '-L/usr/local/lib/other', '-lother')),
        ('-L /usr/local/special -l special',
         ('-L/usr/local/special', '-lspecial'))
    ])
    def test_constructor_link_arguments(self, arg_str: str,
                                        expected: Tuple[str],
                                        user_path: Path, user_pkg_path: Path,
                                        tmp_path: Path):
        """
        Checks linker path arguments are compressed.

        Libraries on user paths.
        """
        (user_pkg_path / 'test.pc').write_text(
            dedent(
                f"""
                Name: test
                Version: 1.0.0
                Description: Many versions.
                URL: http://example.com/version
                Cflags: -Iopt/special/include
                Cflags.private: -I/opt/special/include
                Libs: {arg_str}
                Libs.private:
                """
            )
        )

        test_unit = Package('test')
        assert test_unit.name == 'test'
        assert test_unit.version == (1, 0, 0)
        assert test_unit.compile_arguments == ('-Iopt/special/include',)
        assert test_unit.link_arguments == expected

    @mark.parametrize('link_type, expected', [
        (LinkType.SHARED, '--shared'), (LinkType.STATIC, '--static')
    ])
    def test_constructor_link_type(self, link_type: LinkType, expected: str,
                                   user_path: Path, user_pkg_path: Path,
                                   tmp_path: Path):
        """
        Checks shared and static linking.

        Mock pkg_config acts as though library is on system path.
        """
        (user_pkg_path / 'test.pc').write_text(
            dedent(
                """
                Name: test
                Version: 1.0.0
                Description: Many versions.
                URL: http://example.com/version
                Cflags: -Iopt/special/include
                Cflags.private: -I/opt/special/include/private
                Libs: -Lopt/special/lib -l test
                Libs.private: -L opt/special/lib/private -ltest-private
                """
            )
        )

        test_unit = Package('test', link_type=link_type)
        assert test_unit.name == 'test'
        assert test_unit.version == (1, 0, 0)
        if link_type == LinkType.SHARED:
            assert test_unit.compile_arguments == ('-Iopt/special/include',)
            assert test_unit.link_arguments == ('-Lopt/special/lib', '-ltest')
        elif link_type == LinkType.STATIC:
            assert test_unit.compile_arguments == (
                '-Iopt/special/include', '-I/opt/special/include/private'
            )
            assert test_unit.link_arguments == (
                '-Lopt/special/lib', '-ltest',
                '-Lopt/special/lib/private', '-ltest-private'
            )
        else:
            assert False

    def test_constructor_tree(self, user_path: Path, user_pkg_path: Path,
                              system_path: Path, system_pkg_path: Path,
                              tmp_path: Path):
        """
        Checks recursive dependencies.
        """
        (system_pkg_path / 'system.pc').write_text(
            dedent(
                """
                Name: system
                Version: 1.0.0
                Description: Is a system thing.
                URL: http://example.com/system
                Cflags: -Iusr/local/include
                Libs: -Lusr/local/lib -lsystem
                """
            )
        )

        (user_pkg_path / 'user.pc').write_text(
            dedent(
                """
                Name: user
                Version: 1.0.0
                Description: Depends on system thing.
                URL: http://example.com/user
                Requires: system
                Cflags: -Iopt/user/include
                Libs: -Lopt/user/lib -luser
                """
            )
        )

        test_unit = Package('user')
        assert test_unit.name == 'user'
        assert test_unit.version == (1, 0, 0)
        assert test_unit.compile_arguments == ('-Iopt/user/include',
                                               '-Iusr/local/include')
        assert test_unit.link_arguments == ('-Lopt/user/lib', '-luser',
                                            '-Lusr/local/lib', '-lsystem')

    @mark.parametrize('version, expected',
                      [
                          ('0.1.0', None), ('1.2.2', None),
                          ('1.2.3', (1, 2, 3)), ('1.2.4', (1, 2, 4))
                      ])
    def test_constructor_specification(self, version: str,
                                       expected: Optional[Tuple[int]],
                                       user_pkg_path: Path):
        (user_pkg_path / 'test.pc').write_text(
            dedent(
                f"""
                Name: test
                Version: {version}
                Description: Many versions.
                URL: http://example.com/test
                Cflags: -Iopt/test/include
                Libs: -Lopt/test/lib -ltest
                """
            )
        )
        if expected is None:
            with raises(PackageException):
                _ = Package('test >= 1.2.3')
        else:  # Package should fulfil specification
            test_unit = Package('test >= 1.2.3')
            assert test_unit.name == 'test'
            assert test_unit.version == expected
