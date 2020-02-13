import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed, escapeRegExp, withRegex } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import { User } from '../../Entities/entities';
import SearchWorker from '../../search_worker';

const ListMoleculeRouter = Router();

interface Filters {
  /* Filters */
  version?: string;
  force_fields?: string;
  q: string;
  author: string;
  categories?: string;
  martinize_versions?: string;
  name?: string;
  alias?: string;

  /* Search modificators */
  is_regex?: string;
  combine?: "0" | "false" | "1" | "true";
  skip?: string;
  limit?: string;
}

ListMoleculeRouter.get('/', (req, res) => {
  (async () => {
    const query: nano.MangoQuery = { selector: {}, limit: 25, skip: 0 };

    const { 
      force_fields, 
      version, 
      q: free_text, 
      author, 
      categories, 
      is_regex,
      martinize_versions,
      name: molecule_name,
      alias,
      combine: __as_version_tree__,
      skip,
      limit,
    } = req.query as Partial<Filters>;

    const with_regex = is_regex === "true" || is_regex === "1";
    const bulk_request = __as_version_tree__ === "false" || __as_version_tree__ === "0";

    const selectors: any[] = [];
    if (force_fields) {
      const ff = force_fields.split(',');
      selectors.push({
        force_field: {
          $in: ff
        }
      });
    }
    if (martinize_versions) {
      const vers = martinize_versions.split(',');

      selectors.push({
        martinize_version: {
          $in: vers
        }
      });
    }
    if (version) {
      selectors.push({
        version: withRegex(version, with_regex)
      });
    }
    if (molecule_name) {
      selectors.push({
        name: withRegex(molecule_name, with_regex)
      });
    }
    if (alias) {
      selectors.push({
        alias: withRegex(alias, with_regex)
      });
    }
    if (author) {
      const users_that_match_query = await Database.user.find({
        limit: 99999,
        selector: {
          $or: [
            // @ts-ignore
            { email: withRegex(author, with_regex) },
            { name: withRegex(author, with_regex) },
          ]
        }
      }) as User[];

      const users_id = users_that_match_query.map(u => u.id);

      selectors.push({
        owner: { $in: users_id }
      });
    }
    if (categories) {
      const cat = categories.split(',');
      selectors.push({
        category: {
          $in: cat
        }
      });
    }
    if (free_text) {
      const search_text = with_regex ? author : escapeRegExp(free_text);
      const search_obj = {
        $regex: '(?i)' + search_text
      };

      selectors.push({
        $or: [
          { category: search_obj },
          { comments: search_obj },
          { name: search_obj },
          { alias: search_obj },
          { command_line: search_obj },
          { force_field: search_obj },
        ]
      });
    }
    
    if (selectors.length) {
      query.selector = { $and: selectors };
    }

    if (skip) {
      if (Number(skip) > 0) {
        query.skip = Number(skip);
      }
    }
    if (limit) {
      const l = Number(limit);

      if (l > 0 && l <= 200) {
        query.limit = l;
      }
    }

    const response = await SearchWorker.query(query, bulk_request);
    response.molecules = response.molecules.map(e => sanitize(e));

    res.json(response);
  })().catch(errorCatcher(res));
});

ListMoleculeRouter.all('/', methodNotAllowed('GET'))

export default ListMoleculeRouter;
