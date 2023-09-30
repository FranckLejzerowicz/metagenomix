# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from metagenomix._io_utils import compute_hash


class Soft(object):
    """
    This class is instantiated for very software of the pipeline and
    will contain all the info about it, that can be looked up by finding this
    class after it is registered as value in the `self.softs` dict, which key
    is the software name itself.

    It contains the following info:
        self.config     : copy of all the config
        self.name       : software name
        self.hash       : text to hash (see `get_hash()` in "commands.py")
        self.hashed     : hashed text (see `get_hash()` in "commands.py")
        self.dir        : software output folder
        self.prev       : previously-run software name
        self.scratch    : scratch location where to run this software
        self.params     : default -> user parameters for this software

        # the following attributes are set when the software commands are made
        self.io         : input/output for the movement to scratch
        self.links      :
        self.inputs     :
        self.outputs    :
        self.defaults   :
        self.cmds       :
        self.bash       :
        self.path       :
        self.status     :
        self.tables     :
        self.dirs       :
        self.messages   :
    """

    def __init__(self, config):
        self.config = config
        self.path = []
        self.prev = None
        self.name = ''
        self.hash = ''
        self.hashed = ''
        self.dir = None
        self.scratch = None   # no use of the scratch file system by default
        self.params = dict(config.params)  # init with default params
        self.io = {}
        self.links = {}
        self.inputs = {}
        self.outputs = {}
        self.defaults = {}
        self.cmds = {}
        self.bash = []
        self.status = []
        self.tables = []
        self.dirs = set()
        self.messages = set()

    def set_soft(self, params, path):
        self.prev, self.name = path[-2:]
        self.path = path
        self.params.update(params[self.name])

    def add_status(
            self,
            tech,
            sam_pool,
            dec=None,
            group=None,
            message=None,
            genome=None
    ):
        """A "status" format consists of a list made of 6 elements:
            - `tech`                    : technology
            - `sam_pool`                : sample or pool assembly group
            - `"Done/To do/<tuple>"`    : Done/To do/list
            - `group`                   : co-assembly group (default to None)
            - `message`                 : message to print (default to None)
            - `genome`                  : name of the genome (default to None)
        """
        # 0 or 1 are given for the
        if dec == 0:
            val = 'Done'
        elif dec == 1:
            val = 'To do'
        else:
            if isinstance(dec, list):
                val = tuple(dec)
            elif isinstance(dec, str):
                val = (dec,)
        row = [tech, sam_pool, val]
        row.extend([None if not x else x for x in [group, message, genome]])
        self.status.append(row)

    def get_to_avoid(self) -> set:
        """Define the set of parameters that can be changed without affecting
        the hash value of the software's "after" folder.

        Returns
        -------
        avoid : set
            Names of parameters to be avoided when computing the hash value
        """
        # variables to not account for in the calculation of the hash value
        avoid = {
            'time', 'nodes', 'mem', 'mem_dim', 'env', 'chunks', 'scratch',
            'machine', 'partition', 'cpus', 'skip_samples', 'path', 'binary'}
        # do not account for 'databases' for specific softwares
        if self.name not in ['filtering', 'databases']:
            avoid.add('databases')
        # do not account for (search) 'terms' for search_* softwares
        if 'search_' in self.name:
            avoid.add('terms')
        return avoid

    def get_hash(self, params, softs):
        """Get info to hash and hash it.
            - `self.soft.hash`      : info to hash
            - `self.soft.hashed`    : hashed info
        """
        avoid = self.get_to_avoid()
        # get all info to hash
        params_dict = dict(x for x in params.items() if x[0] not in avoid)
        hashes = []
        for r in range(2, len(self.path)):
            hashes.append(softs[tuple(self.path[1:r])])
        self.hash = (params_dict, self.path, hashes)
        # hash all this info
        self.hashed = compute_hash(self.hash)
        # fill lookup dict with  hash of current path to re-use for next path
        softs[tuple(self.path[1:])] = (params_dict, self.path)

    # def add_to_path(self, softs):
    #     # if self.prev is None:
    #     print()
    #     print(self.prev)
    #     print(self.name)
    #     print(softs)
    #     if self.prev == 'None':
    #         self.path = ['fastq', self.name]
    #     else:
    #         self.path = softs[(self.prev)].path + [self.name]
    #     print("self.path:", self.path)

