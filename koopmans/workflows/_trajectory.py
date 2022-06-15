


from ._workflow import Workflow
from typing import List
from ase import Atoms, io
import json as json_ext
from koopmans import utils


load_results_from_output = True


class TrajectoryWorkflow(Workflow):

    
    def __init__(self, *args, **kwargs):
        self.snapshots: List[Atoms] = kwargs.pop('snapshots', [])
        super().__init__(*args, **kwargs)


    def _run(self) -> None:
         # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import KoopmansDFPTWorkflow, KoopmansDSCFWorkflow

        self.n_max = 4
        
        for snapshot in self.snapshots:
        
            self.atoms.set_positions(snapshot.positions)
            if self.parameters.method == 'dfpt':
                workflow = KoopmansDFPTWorkflow(**self.wf_kwargs)
                self.run_subworkflow(workflow)
            else:
                dscf_workflow = KoopmansDSCFWorkflow(**self.wf_kwargs)
                self.run_subworkflow(dscf_workflow)

    @classmethod
    def _fromjsondct(cls, bigdct: str):
        # print("INside of _from_json")
        # assert(False)

        snapshots_file = bigdct['setup'].pop('snapshots')
        snapshots = io.read(snapshots_file) #TODO raise ValueError if not provided
        if isinstance(snapshots, Atoms):
            snapshots = [snapshots]
        bigdct['setup']['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])

        wf = super(TrajectoryWorkflow, cls)._fromjsondct(bigdct)

        wf.snapshots = snapshots
        return wf