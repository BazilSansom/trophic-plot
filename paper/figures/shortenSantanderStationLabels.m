function labelsOut = shortenSantanderStationLabels(labelsIn)
%SHORTENSANTANDERSTATIONLABELS_  Reader-friendly display labels for TfL stations.
% Display-only relabelling: does NOT merge nodes.

labelsOut = cellstr(string(labelsIn));

from = { ...
    'Argyle Street, Kings Cross';
    'Holborn Circus, Holborn';
    'Newgate Street, St. Paul''s';
    'Newgate Street , St. Paul''s';
    'Queen Street 1, Bank';
    'Queen Street 2, Bank';
    'Brushfield Street, Liverpool Street';
    'Finsbury Circus, Liverpool Street';
    'Wormwood Street, Liverpool Street';
    'St. James''s Square, St. James''s';
    'Smith Square, Westminster';
    'Baylis Road, Waterloo';
    'Waterloo Station 1, Waterloo';
    'Waterloo Station 2, Waterloo';
    'Waterloo Station 3, Waterloo';
    'Eastbourne Mews, Paddington';
    'London Street, Paddington';
    'Soho Square, Soho';
    'Soho Square , Soho';
    'Moorfields, Moorgate';
    'Battersea Park Road, Nine Elms';
    'Duke Street Hill, London Bridge';
    'Crosswall, Tower';
    'Cheapside, Bank';
    'Old Street Station, St. Luke''s';
    'Tooley Street, Bermondsey';
    'Hop Exchange, The Borough';
    'Little Argyll Street, West End';
    'Hyde Park Corner, Hyde Park';
    'Royal London Hospital, Whitechapel';
    'Blackfriars Station, St. Paul''s';
    'Warren Street Station, Euston';
    'Howick Place, Westminster';
    'William IV Street, Strand';
    'Black Lion Gate, Kensington Gardens';
    'Exhibition Road, Knightsbridge';
    'Allington Street, Victoria';
    'Warwick Row, Westminster';
    'Malet Street, Bloomsbury';
    'Great Tower Street, Monument';
    'Albert Gate, Hyde Park';
    'Westminster Pier, Westminster';
    'London Fields, Hackney Central';
    'Derry Street, Kensington'
    'Strata, Elephant & Castle' 
    'Worship Street, Shoreditch'
    'Lambeth Palace Road, Waterloo'
    'Imperial College, Knightsbridge'
    'Howland Street, Fitzrovia'
    'Boston Place, Marylebone'
    'Wellington Street , Strand'
    'Braham Street, Aldgate'
    'Eccleston Place, Victoria'
    'Bunhill Row, Moorgate'
    'Drayton Gardens, West Chelsea'
    'Museum of London, Barbican'
    'Craven Street, Strand'
    };

to = { ...
    'Kings Cross';
    'Holborn';
    'St Paul''s';
    'St Paul''s';
    'Bank 1';
    'Bank 2';
    'Liverpool St (BS)';
    'Liverpool St (FC)';
    'Liverpool St (WS)';
    'St James''s';
    'Westminster (Smith Sq)';
    'Waterloo (Baylis)';
    'Waterloo 1';
    'Waterloo 2';
    'Waterloo 3';
    'Paddington (Eastbourne)';
    'Paddington (London St)';
    'Soho';
    'Soho';
    'Moorgate';
    'Nine Elms';
    'London Bridge';
    'Tower';
    'Bank 3';
    'Old Street';
    'Bermondsey';
    'The Borough';
    'Little Argyll Street, West End';
    'Hyde Park Corner';
    'Whitechapel';
    'Blackfriars';
    'Warren Street';
    'Westminster (H.Place)';
    'Strand';
    'Kensington Gardens';
    'Knightsbridge';
    'Victoria';
    'Westminster (Warwick Row)';
    'Bloomsbury';
    'Monument';
    'Hyde Park (Albert Gate)';
    'Westminster Pier';
    'Hackney Central';
    'Kensington (Derry St)';
    'Elephant & Castle';
    'Shoreditch';
    'Waterloo (LPlce Rd)';
    'Knightsbridge (Imperial)';
    'Fitzrovia (Howland St)';
    'Marylebone (Boston Pl)';
    'Strand (Wellington St)';
    'Aldgate (Braham St)';
    'Victoria (E.Pl)';
    'Moorgate (B.Row)';
    'West Chelsea';
    'Barbican';
    'Strand (Crvn.St)';
    };

assert(numel(from) == numel(to), 'Label map lengths differ.');

for i = 1:numel(from)
    hit = strcmpi(labelsOut, from{i});
    labelsOut(hit) = to(i);
end
end