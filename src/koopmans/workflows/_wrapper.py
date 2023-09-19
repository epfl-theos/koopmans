from typing import Any, Dict, Generic, List, Type, TypeVar

from ._workflow import Workflow

WrapType = TypeVar('WrapType', bound='WorkflowWrapper')
WorkflowType = TypeVar('WorkflowType', bound='Workflow')


class WorkflowWrapper(Workflow):
    _wrapper_args: List[str] = []
    subworkflow_classes: List[Workflow]

    def __init__(self, subworkflow_classes: List[Workflow], **kwargs):
        self.subworkflow_classes = subworkflow_classes
        for arg in self._wrapper_args:
            setattr(self, arg, kwargs.pop(arg, None))
        super().__init__(**kwargs)

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
        for kwarg in cls._wrapper_args + ['_subworkflow_class']:
            if kwarg in dct:
                kwargs[kwarg] = dct.pop(kwarg.lstrip('_'))
        return super(WorkflowWrapper, cls).fromdict(dct, **kwargs)

    @classmethod
    def fromparent(cls: Type[WorkflowType], parent: Workflow, **kwargs) -> WorkflowType:
        if isinstance(parent, WorkflowWrapper):
            kwargs['subworkflow_classes'] = parent.subworkflow_classes[:-1]
        return super(WorkflowWrapper, cls).fromparent(parent, **kwargs)

    @property
    def subworkflow_class(self) -> Workflow:
        return self.subworkflow_classes[-1]


def WorkflowWrapperFactory(wrapper: Generic[WrapType], subworkflow: Workflow, **kwargs) -> Generic[WrapType]:
    # Function that wraps subworkflow with a workflow wrapper, co-opting the fromparent method to use the subworkflow
    # settings to initialize the wrapped workflow

    # Initialize without subworkflow_classes because these will be added by the wrapper
    wf = wrapper.fromparent(subworkflow, subworkflow_classes=[], **kwargs)

    # fromparent, by default, sets wf.parent to subworkflow whereas we want to set it to None
    wf.parent = None

    # Adding the subworkflow_classes
    if isinstance(subworkflow, WorkflowWrapper):
        wf.subworkflow_classes = subworkflow.subworkflow_classes + [subworkflow.__class__]
    else:
        wf.subworkflow_classes = [subworkflow.__class__]

    return wf
