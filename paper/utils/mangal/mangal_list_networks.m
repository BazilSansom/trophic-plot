function nets = mangal_list_networks(dataset_id)
%MANGAL_LIST_NETWORKS  Return *all* networks for a dataset_id (paged).

limit = 200;
offset = 0;
nets = struct([]);

while true
    page = mangal_api_get("network", "dataset_id", dataset_id, "limit", limit, "offset", offset);
    if isempty(page)
        break;
    end
    nets = [nets; page(:)]; %#ok<AGROW>
    offset = offset + limit;
end
end