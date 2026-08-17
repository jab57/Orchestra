"""
Safeguard for Orchestra experiment scripts: verifies a tool response is actually
about the gene that was requested before any downstream code trusts it.

Built after a pilot investigation found causal_chain_analysis's effector/scaffold
routing (PPI -> TF-partner substitution, e.g. CCND1 -> ESR1) could silently return
a different gene's data than what was requested. That specific routing does not
exist in _run_network_comparison_path (compare_network_contexts) -- verified by
reading the code -- but this guard is cheap defense-in-depth against any tool
whose behavior isn't as thoroughly verified.
"""


class GeneMismatchError(Exception):
    """Raised when a tool's response payload claims to be about a different gene
    than what was requested -- data integrity failure, must not be silently used."""


def assert_gene_match(result: dict, requested_gene: str, *, field: str = "gene") -> None:
    """
    Verify a RegNetAgents/CASCADE tool response is actually about the gene we asked for.

    Call this immediately after every tool call in the data-collection pipeline, before
    the result is used for anything (candidate extraction, scoring, storage). Raises
    rather than warns -- a silently-accepted wrong-gene result corrupts downstream
    statistics in a way that's very hard to detect after the fact.
    """
    returned_gene = result.get(field)
    if returned_gene is None:
        raise GeneMismatchError(
            f"Response has no top-level {field!r} field to verify against "
            f"requested_gene={requested_gene!r}. Response keys: {list(result.keys())}"
        )
    if returned_gene.upper() != requested_gene.upper():
        raise GeneMismatchError(
            f"Gene mismatch: requested {requested_gene!r} but response is about "
            f"{returned_gene!r}. Do not use this result -- discard and re-query, "
            f"or if this is expected (e.g. effector/TF-partner routing), handle explicitly "
            f"rather than silently accepting it as the requested gene's data."
        )
