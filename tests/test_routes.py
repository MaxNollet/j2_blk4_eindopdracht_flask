def test_home_route_status_code(client):
    """A function which tests the home-route."""
    route = "/"
    rv = client.get(route)
    assert rv.status_code == 200


def test_query_route_status_code(client):
    """A function which tests the query-route."""
    route = "/query"
    rv = client.get(route)
    assert rv.status_code == 200
