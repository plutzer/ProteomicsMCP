"""
Test script for PhosphoSitePlus MCP tools
"""
import logging
import sys

# Import the tools directly
from psp_proteomics import (
    get_kinase_substrates,
    get_regulatory_sites,
    get_disease_sites,
    find_upstream_kinases,
    get_interaction_network,
    get_pathway_context,
    get_evidence_summary
)

def print_divider(title):
    print('\n' + '='*70)
    print(f'  {title}')
    print('='*70)

def print_csv_preview(csv_string, max_rows=5):
    """Print first few rows of a CSV string."""
    lines = csv_string.split('\n')
    for i, line in enumerate(lines[:max_rows]):
        print(f"  {line}")
    if len(lines) > max_rows:
        print(f"  ... ({len(lines) - max_rows} more rows)")

def test_kinase_substrates():
    print_divider("Test 1: get_kinase_substrates")
    print("Query: AKT1,MTOR,GSK3B")

    try:
        result = get_kinase_substrates('AKT1,MTOR,GSK3B')

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        print(f"\nFound data for {len(result)} kinases")
        for kinase, csv_data in result.items():
            print(f"\n{kinase}:")
            print_csv_preview(csv_data, 5)

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_regulatory_sites():
    print_divider("Test 2: get_regulatory_sites")
    print("Query: TP53,AKT1_S473")

    try:
        result = get_regulatory_sites('TP53,AKT1_S473')

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        print(f"\nFound data for {len(result)} queries")
        for query_item, csv_data in result.items():
            print(f"\n{query_item}:")
            print_csv_preview(csv_data, 5)

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_disease_sites():
    print_divider("Test 3: get_disease_sites")
    print("Query: CTNNB1,TP53")

    try:
        result = get_disease_sites('CTNNB1,TP53')

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        print(f"\nFound data for {len(result)} queries")
        for query_item, csv_data in result.items():
            print(f"\n{query_item}:")
            print_csv_preview(csv_data, 5)

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_upstream_kinases():
    print_divider("Test 4: find_upstream_kinases")
    print("Query: AKT1_S473,GSK3B_S9,MAPK1")

    try:
        result = find_upstream_kinases('AKT1_S473,GSK3B_S9,MAPK1')

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        print(f"\nFound data for {len(result)} queries")
        for query_item, csv_data in result.items():
            print(f"\n{query_item}:")
            print_csv_preview(csv_data, 5)

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_interaction_network():
    print_divider("Test 5: get_interaction_network")
    print("Query: MTOR_S2448,AKT1 (depth=1)")

    try:
        result = get_interaction_network('MTOR_S2448,AKT1', depth=1)

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        print(f"\nNetwork summary: {result['summary']}")
        print("\nNodes:")
        print_csv_preview(result['nodes'], 10)
        print("\nEdges:")
        print_csv_preview(result['edges'], 10)

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_pathway_context():
    print_divider("Test 6: get_pathway_context")
    print("Query: AKT1_S473")

    try:
        result = get_pathway_context('AKT1_S473')

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        print(f"\nQuery: {result['query']}")
        print(f"Summary: {result['summary']}")

        if 'upstream' in result:
            print("\nUpstream kinases:")
            print_csv_preview(result['upstream'], 5)

        if 'downstream' in result:
            print("\nDownstream substrates:")
            print_csv_preview(result['downstream'], 5)

        if 'functional_outcomes' in result:
            print("\nFunctional outcomes:")
            print_csv_preview(result['functional_outcomes'], 5)

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_evidence_summary():
    print_divider("Test 7: get_evidence_summary")
    print("Query: MTOR_S2448")

    try:
        result = get_evidence_summary('MTOR_S2448')

        if 'error' in result:
            print(f"ERROR: {result['error']}")
            return False

        for query_item, data in result.items():
            print(f"\nQuery: {query_item}")
            print(f"Summary: {data['summary']}")

            if 'top_kinases' in data:
                print("\nTop kinases:")
                print_csv_preview(data['top_kinases'], 5)

            if 'functional_categories' in data:
                print("\nFunctional categories:")
                print_csv_preview(data['functional_categories'], 5)

            if 'disease_associations' in data:
                print("\nDisease associations:")
                print_csv_preview(data['disease_associations'], 5)

            if 'key_pmids' in data:
                print(f"\nKey PMIDs: {data['key_pmids'][:100]}...")

        return True
    except Exception as e:
        print(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    print("\n" + "="*70)
    print("  PhosphoSitePlus MCP Tools Test Suite")
    print("="*70)

    tests = [
        ("Kinase Substrates", test_kinase_substrates),
        ("Regulatory Sites", test_regulatory_sites),
        ("Disease Sites", test_disease_sites),
        ("Upstream Kinases", test_upstream_kinases),
        ("Interaction Network", test_interaction_network),
        ("Pathway Context", test_pathway_context),
        ("Evidence Summary", test_evidence_summary),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            print(f"\nFATAL ERROR in {test_name}: {e}")
            results.append((test_name, False))

    # Summary
    print_divider("Test Results Summary")
    passed_count = 0
    for test_name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {test_name:30s}: {status}")
        if passed:
            passed_count += 1

    print(f"\nTotal: {passed_count}/{len(tests)} tests passed")

    if passed_count == len(tests):
        print("\nAll tests passed! Ready for MCP integration.")
        return 0
    else:
        print("\nSome tests failed. Please review errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
