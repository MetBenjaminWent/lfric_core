##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Makes use of pkg-config files to find out about libraries.
"""
from os import environ as os_environ
from enum import StrEnum
from re import match as re_match
from subprocess import run, PIPE
from typing import Dict, Iterable, List, Tuple, Union


class LinkType(StrEnum):
    SHARED = '--shared'
    STATIC = '--static'


class PackageException(Exception):
    pass


class Package:
    """
    Holds details of a library.
    """
    def __init__(self, specification: str,
                 link_type: LinkType = LinkType.SHARED):
        """
        Constructs Package object from details held in pkg-config files.

        Version requirement may be specified in the same way pkg-config
        understands them. i.e. =, >, <, >= and <=

        :param specification: Package name with optional version requirement.
        :raises PackageException: A package fulfilling the specification was
                not found.
        """
        match = re_match(r'(\w+)[ \t<>=]?', specification)
        if match:
            self.__name = match.group(1)
        else:
            raise PackageException(
                "Unable to parse specification: " + specification
            )

        self.__pkg_config(specification, ['--print-errors', '--exists'])

        version = self.__pkg_config(specification, ['--modversion'])
        if version:
            self.__version = tuple([int(component)
                                    if component.isdigit() else component
                                    for component in version[0].split('.')])
        else:  # No version string
            self.__version = tuple()

        compile_arguments = self.__pkg_config(specification,
                                              ['--cflags', link_type])
        self.__compile_arguments = self.__split_arguments(compile_arguments)

        # ToDo: It may be necessary to use --pure with static link type.
        #
        link_arguments = self.__pkg_config(specification,
                                           ['--libs', link_type])
        self.__link_arguments = self.__split_arguments(link_arguments)

    @staticmethod
    def __pkg_config(specification: str, arguments: List[str]) -> List[str]:
        environment: Dict[str, str] = {'PKG_CONFIG_DEBUG_SPEW': 'YES'}
        if 'PKG_CONFIG_LIBDIR' in os_environ:
            environment['PKG_CONFIG_LIBDIR'] = os_environ['PKG_CONFIG_LIBDIR']
        if 'PKG_CONFIG_PATH' in os_environ:
            environment['PKG_CONFIG_PATH'] = os_environ['PKG_CONFIG_PATH']
        command = ['pkg-config']
        command.extend(arguments)
        command.append(specification)
        process = run(command, stdout=PIPE, stderr=PIPE, encoding='utf8',
                      env=environment)
        if process.returncode != 0:
            raise PackageException(
                f"Failed to run [{command}]: {process.stderr}"
            )
        if process.stderr:
            raise PackageException(process.stderr)
        return process.stdout.split()

    @staticmethod
    def __split_arguments(arguments: Iterable[str]) -> Tuple[str, ...]:
        result: List[str] = []
        directive = ''
        for argument in arguments:
            if argument in ['-I', '-L', '-l']:
                directive = argument
                continue
            if directive:
                result.append(directive + argument)
                directive = ''
                continue
            result.append(argument)
        return tuple(result)

    @property
    def name(self) -> str:
        """
        Gets library name.
        """
        return self.__name

    @property
    def version(self) -> Tuple[Union[int, str], ...]:
        """
        Gets library version.
        """
        return self.__version

    @property
    def compile_arguments(self) -> Tuple[str, ...]:
        """
        Gets arguments for compiler.

        Arguments are canonicalised into "no space" form.
        """
        return self.__compile_arguments

    @property
    def link_arguments(self) -> Tuple[str, ...]:
        """
        Gets arguments for linker.

        Arguments are canonicalised into "no space" form.
        """
        return self.__link_arguments
